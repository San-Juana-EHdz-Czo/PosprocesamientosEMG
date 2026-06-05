# =========================================================
# LIBRERÍAS
# =========================================================

library(biosignalEMG)
library(dplyr)
library(knitr)
library(kableExtra)

# =========================================================
# CONTADOR DE TIEMPO
# =========================================================

tiempo_inicio <- Sys.time()

# =========================================================
# PARÁMETROS GENERALES
# =========================================================

n_sim <- 1000

fs <- 1000

# =========================================================
# PARÁMETROS ÓPTIMOS
# =========================================================

# Operadores morfológicos
m <- 5

# Duración mínima (~30 ms)
min_len <- 30

# Majority filter
majority_window <- 5

# Ising
lambda_ising <- 1.5
gamma_ising  <- 2

# Número de iteraciones
max_iter_greedy <- 50
n_iter_metro    <- 150
n_iter_gibbs    <- 100

# =========================================================
# FUNCIONES AUXILIARES
# =========================================================

# ---------------------------------------------------------
# RUNS
# ---------------------------------------------------------

count_runs <- function(binary_signal) {
  1 + sum(diff(binary_signal) != 0)
}

compute_Z <- function(env, t) {
  
  n <- length(env)
  
  b_hat <- ifelse(env >= t, 1, 0)
  
  R_t <- count_runs(b_hat)
  
  p_t <- mean(env >= t)
  q_t <- 1 - p_t
  
  if (p_t == 0 || p_t == 1) return(Inf)
  
  mu_R <- 1 + 2 * (n - 1) * p_t * q_t
  
  var_R <- 2 * p_t * q_t *
    (2 * n - 3 - 2 * (3 * n - 5) * p_t * q_t)
  
  if (var_R <= 0) return(Inf)
  
  (R_t - mu_R) / sqrt(var_R)
}

find_threshold_runs <- function(env,
                                num_thresholds = 50) {
  
  t_values <- seq(
    min(env),
    max(env),
    length.out = num_thresholds
  )
  
  Z_values <- sapply(
    t_values,
    function(t) compute_Z(env, t)
  )
  
  best_index <- which.min(Z_values)
  
  list(
    threshold = t_values[best_index]
  )
}

runs_threshold_method <- function(signal) {
  
  env <- signal$values
  
  res <- find_threshold_runs(env)
  
  t_opt <- res$threshold
  
  ifelse(env >= t_opt, 1, 0)
}

# ---------------------------------------------------------
# BAYES
# ---------------------------------------------------------

em_gaussian_mixture <- function(x,
                                max_iter = 100,
                                tol = 1e-5) {
  
  sigma0 <- sd(x) / 2
  sigma1 <- sd(x)
  
  lambda0 <- 0.5
  lambda1 <- 0.5
  
  loglik_old <- -Inf
  
  for (iter in 1:max_iter) {
    
    f0 <- dnorm(x, mean = 0, sd = sigma0)
    f1 <- dnorm(x, mean = 0, sd = sigma1)
    
    w0 <- lambda0 * f0 /
      (lambda0 * f0 + lambda1 * f1)
    
    w1 <- 1 - w0
    
    lambda0 <- mean(w0)
    lambda1 <- 1 - lambda0
    
    sigma0 <- sqrt(sum(w0 * x^2) / sum(w0))
    sigma1 <- sqrt(sum(w1 * x^2) / sum(w1))
    
    loglik <- sum(log(lambda0 * f0 + lambda1 * f1))
    
    if (abs(loglik - loglik_old) < tol) break
    
    loglik_old <- loglik
  }
  
  list(
    sigma0 = sigma0,
    sigma1 = sigma1
  )
}

compute_bayes_threshold <- function(params,
                                    gamma = 1) {
  
  sigma0 <- params$sigma0
  sigma1 <- params$sigma1
  
  if (sigma0 == sigma1)
    return(mean(c(sigma0, sigma1)))
  
  C <- gamma * (sigma1 / sigma0)
  
  sqrt(2 * log(C)) *
    (sigma0 * sigma1) /
    abs(sigma1 - sigma0)
}

bayesian_threshold_method <- function(signal,
                                      gamma = 1) {
  
  env <- signal$values
  
  params <- em_gaussian_mixture(env)
  
  t <- compute_bayes_threshold(
    params,
    gamma
  )
  
  ifelse(env >= t, 1, 0)
}

# =========================================================
# POSPROCESAMIENTO
# =========================================================

# ---------------------------------------------------------
# DILATACIÓN
# ---------------------------------------------------------

dilation <- function(b, m) {
  
  n <- length(b)
  out <- numeric(n)
  
  for (i in 1:n) {
    
    idx <- max(1, i - m):min(n, i + m)
    
    out[i] <- max(b[idx])
  }
  
  out
}

# ---------------------------------------------------------
# EROSIÓN
# ---------------------------------------------------------

erosion <- function(b, m) {
  
  n <- length(b)
  out <- numeric(n)
  
  for (i in 1:n) {
    
    idx <- max(1, i - m):min(n, i + m)
    
    out[i] <- min(b[idx])
  }
  
  out
}

# ---------------------------------------------------------
# OPENING / CLOSING
# ---------------------------------------------------------

opening <- function(b, m) {
  dilation(erosion(b, m), m)
}

closing <- function(b, m) {
  erosion(dilation(b, m), m)
}

# ---------------------------------------------------------
# REMOVE SHORT ACTIVITY
# ---------------------------------------------------------

remove_short_activity <- function(b,
                                  min_len) {
  
  r <- rle(b)
  
  idx <- which(
    r$values == 1 &
      r$lengths < min_len
  )
  
  r$values[idx] <- 0
  
  inverse.rle(r)
}

# ---------------------------------------------------------
# REMOVE SHORT SILENCE
# ---------------------------------------------------------

remove_short_silence <- function(b,
                                 min_len) {
  
  r <- rle(b)
  
  idx <- which(
    r$values == 0 &
      r$lengths < min_len
  )
  
  r$values[idx] <- 1
  
  inverse.rle(r)
}

# ---------------------------------------------------------
# PIPELINE COMPLETO
# ---------------------------------------------------------

postprocess <- function(b,
                        m,
                        min_len) {
  
  b_close <- closing(b, m)
  
  b_open <- opening(b_close, m)
  
  b_act <- remove_short_activity(
    b_open,
    min_len
  )
  
  remove_short_silence(
    b_act,
    min_len
  )
}

# =========================================================
# MODELOS DE ISING
# =========================================================

# ---------------------------------------------------------
# GREEDY
# ---------------------------------------------------------

ising_denoise <- function(y,
                          lambda = 1,
                          gamma = 1,
                          max_iter = 50) {
  
  x <- y
  n <- length(y)
  
  for (iter in 1:max_iter) {
    
    changed <- FALSE
    
    for (i in 1:n) {
      
      left  <- ifelse(i > 1, x[i - 1], 0)
      right <- ifelse(i < n, x[i + 1], 0)
      
      e0 <- lambda *
        (abs(0 - left) + abs(0 - right)) +
        gamma * abs(0 - y[i])
      
      e1 <- lambda *
        (abs(1 - left) + abs(1 - right)) +
        gamma * abs(1 - y[i])
      
      new_val <- ifelse(e1 < e0, 1, 0)
      
      if (new_val != x[i]) {
        x[i] <- new_val
        changed <- TRUE
      }
    }
    
    if (!changed) break
  }
  
  x
}

# ---------------------------------------------------------
# METROPOLIS
# ---------------------------------------------------------

ising_metropolis <- function(y,
                             lambda = 1,
                             gamma = 1,
                             n_iter = 1000) {
  
  x <- y
  n <- length(y)
  
  energy <- function(z) {
    
    lambda * sum(abs(diff(z))) +
      gamma * sum(abs(z - y))
  }
  
  E_current <- energy(x)
  
  for (iter in 1:n_iter) {
    
    i <- sample.int(n, 1)
    
    x_new <- x
    x_new[i] <- 1 - x_new[i]
    
    E_new <- energy(x_new)
    
    dE <- E_new - E_current
    
    if (dE <= 0 ||
        runif(1) < exp(-dE)) {
      
      x <- x_new
      E_current <- E_new
    }
  }
  
  x
}

# ---------------------------------------------------------
# GIBBS
# ---------------------------------------------------------

ising_gibbs <- function(y,
                        lambda = 1,
                        gamma = 1,
                        n_iter = 500) {
  
  x <- y
  n <- length(y)
  
  for (iter in 1:n_iter) {
    
    for (i in 1:n) {
      
      left  <- ifelse(i > 1, x[i - 1], 0)
      right <- ifelse(i < n, x[i + 1], 0)
      
      e1 <- lambda *
        (abs(1 - left) + abs(1 - right)) +
        gamma * abs(1 - y[i])
      
      e0 <- lambda *
        (abs(0 - left) + abs(0 - right)) +
        gamma * abs(0 - y[i])
      
      p1 <- exp(-e1)
      p0 <- exp(-e0)
      
      prob <- p1 / (p1 + p0)
      
      x[i] <- rbinom(1, 1, prob)
    }
  }
  
  x
}

# =========================================================
# MAJORITY FILTER
# =========================================================

majority_filter <- function(y,
                            window = 5) {
  
  stats::filter(
    y,
    rep(1 / window, window),
    sides = 2
  ) |> 
    round() |>
    as.numeric() |>
    (\(x){
      x[is.na(x)] <- y[is.na(x)]
      x
    })()
}

# =========================================================
# PROGRAMACIÓN DINÁMICA
# =========================================================

binary_segmentation_dp <- function(y,
                                   lambda = 1){
  
  n <- length(y)
  
  cost <- rep(Inf, n + 1)
  
  prev <- rep(0, n + 1)
  
  cost[1] <- 0
  
  for(i in 2:(n + 1)){
    
    for(j in 1:(i - 1)){
      
      segment <- y[j:(i - 1)]
      
      val <- round(mean(segment))
      
      seg_cost <- sum(abs(segment - val))
      
      total_cost <- cost[j] +
        seg_cost +
        lambda
      
      if(total_cost < cost[i]){
        
        cost[i] <- total_cost
        
        prev[i] <- j
      }
    }
  }
  
  x <- rep(0, n)
  
  i <- n + 1
  
  while(i > 1){
    
    j <- prev[i]
    
    segment <- y[j:(i - 1)]
    
    val <- round(mean(segment))
    
    x[j:(i - 1)] <- val
    
    i <- j
  }
  
  x
}

# =========================================================
# REGULARIZACIÓN PROBABILÍSTICA
# =========================================================

probabilistic_regularization <- function(y,
                                         alpha = 0.7) {
  
  n <- length(y)
  
  x <- y
  
  for(i in 2:(n - 1)){
    
    neighborhood <- mean(
      c(y[i - 1], y[i], y[i + 1])
    )
    
    x[i] <- ifelse(
      neighborhood >= alpha,
      1,
      0
    )
  }
  
  x
}

# =========================================================
# CONTENEDOR
# =========================================================

n_classifiers <- 4
n_methods <- 13

results <- vector(
  "list",
  n_sim * n_classifiers * n_methods
)

counter <- 1

# =========================================================
# SIMULACIONES
# =========================================================

for(sim in 1:n_sim){
  
  cat(
    "Simulación:",
    sim,
    "de",
    n_sim,
    "\n"
  )
  
  # -------------------------------------------------------
  # SEÑAL
  # -------------------------------------------------------
  
  emgx <- syntheticemg(
    n.length.out = 5000,
    on.sd = 1,
    on.duration.mean = 350,
    on.duration.sd = 10,
    off.sd = 0.05,
    off.duration.mean = 300,
    off.duration.sd = 20,
    on.mode.pos = 0.75,
    shape.factor = 0.5,
    samplingrate = fs,
    units = "mV",
    data.name = "Synthetic EMG"
  )
  
  truth <- emgx$on.off
  
  # -------------------------------------------------------
  # PREPROCESAMIENTO
  # -------------------------------------------------------
  
  emgx_rect <- rectification(
    emgx,
    rtype = "fullwave"
  )
  
  env_emgx <- envelope(
    emgx_rect,
    method = "RMS",
    wsize = 50
  )
  
  # -------------------------------------------------------
  # CLASIFICADORES
  # -------------------------------------------------------
  
  classifiers <- list(
    
    bonato =
      onoff_bonato(
        env_emgx,
        sigma_n = 1
      ),
    
    singlethres =
      onoff_singlethres(
        env_emgx
      ),
    
    threshold_runs =
      runs_threshold_method(
        env_emgx
      ),
    
    threshold_bayes =
      bayesian_threshold_method(
        env_emgx
      )
  )
  
  # -------------------------------------------------------
  # POSPROCESAMIENTO
  # -------------------------------------------------------
  
  for(clf_name in names(classifiers)){
    
    base_signal <- classifiers[[clf_name]]
    
    if(is.list(base_signal))
      base_signal <- base_signal$binary
    
    post_methods <- list(
      
      original =
        base_signal,
      
      dilation =
        dilation(base_signal, m),
      
      erosion =
        erosion(base_signal, m),
      
      closing =
        closing(base_signal, m),
      
      opening =
        opening(base_signal, m),
      
      remove_short_activity =
        remove_short_activity(
          base_signal,
          min_len
        ),
      
      remove_short_silence =
        remove_short_silence(
          base_signal,
          min_len
        ),
      
      all =
        postprocess(
          base_signal,
          m,
          min_len
        ),
      
      ising_greedy =
        ising_denoise(
          base_signal,
          lambda = lambda_ising,
          gamma = gamma_ising,
          max_iter = max_iter_greedy
        ),
      
      ising_metropolis =
        ising_metropolis(
          base_signal,
          lambda = lambda_ising,
          gamma = gamma_ising,
          n_iter = n_iter_metro
        ),
      
      ising_gibbs =
        ising_gibbs(
          base_signal,
          lambda = lambda_ising,
          gamma = gamma_ising,
          n_iter = n_iter_gibbs
        ),
      
      majority_filter =
        majority_filter(
          base_signal,
          window = majority_window
        ),
      
      probabilistic_regularization =
        probabilistic_regularization(
          base_signal,
          alpha = 0.7
        )
    )
    
    # -----------------------------------------------------
    # MÉTRICAS
    # -----------------------------------------------------
    
    for(method_name in names(post_methods)){
      
      pred <- post_methods[[method_name]]
      
      andp_val <- ANDP(pred, truth)
      
      mnchpd_val <- MNChPD(pred, truth)
      
      pce_val <- PCE(pred, truth)
      
      # PR devuelve TPR y FPR
      pr_val <- PR(pred, truth, t=50)
      
      # Algunas versiones lo regresan como vector nombrado
      tpr_val <- pr_val["TPR"]
      
      fpr_val <- pr_val["FPR"]
      
      results[[counter]] <- data.frame(
        
        simulation = sim,
        
        classifier = clf_name,
        
        postprocess = method_name,
        
        ANDP = andp_val,
        
        MNChPD = mnchpd_val,
        
        PCE = pce_val,
        
        TPR = tpr_val,
        
        FPR = fpr_val
      )
      
      counter <- counter + 1
    }
  }
}

# =========================================================
# DATAFRAME FINAL
# =========================================================

results_df <- bind_rows(results)

stopifnot(
  !any(is.na(results_df$postprocess))
)

# =========================================================
# RESUMEN ESTADÍSTICO
# =========================================================

summary_results <- results_df %>%
  
  group_by(
    classifier,
    postprocess
  ) %>%
  
  summarise(
    
    ANDP_mean = mean(
      ANDP,
      na.rm = TRUE
    ),
    
    ANDP_sd = sd(
      ANDP,
      na.rm = TRUE
    ),
    
    MNChPD_mean = mean(
      MNChPD,
      na.rm = TRUE
    ),
    
    MNChPD_sd = sd(
      MNChPD,
      na.rm = TRUE
    ),
    
    PCE_mean = mean(
      PCE,
      na.rm = TRUE
    ),
    
    PCE_sd = sd(
      PCE,
      na.rm = TRUE
    ),
    
    TPR_mean = mean(
      TPR,
      na.rm = TRUE
    ),
    
    TPR_sd = sd(
      TPR,
      na.rm = TRUE
    ),
    
    FPR_mean = mean(
      FPR,
      na.rm = TRUE
    ),
    
    FPR_sd = sd(
      FPR,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

# =========================================================
# TABLA BONITA
# =========================================================

tabla_final <- summary_results %>%
  
  mutate(
    
    ANDP = sprintf(
      "%.2f ± %.2f",
      ANDP_mean,
      ANDP_sd
    ),
    
    MNChPD = sprintf(
      "%.2f ± %.2f",
      MNChPD_mean,
      MNChPD_sd
    ),
    
    PCE = sprintf(
      "%.2f ± %.2f",
      PCE_mean,
      PCE_sd
    ),
    TPR = sprintf(
      "%.2f ± %.2f",
      TPR_mean,
      TPR_sd
    ),
    
    FPR = sprintf(
      "%.2f ± %.2f",
      FPR_mean,
      FPR_sd
    )
  ) %>%
  
  select(
    classifier,
    postprocess,
    ANDP,
    MNChPD,
    PCE,
    TPR,
    FPR
  )

kable(
  tabla_final,
  
  caption = "Resultados promedio de las métricas para cada clasificador y método de posprocesamiento",
  
  col.names = c(
    "Clasificador",
    "Posprocesamiento",
    "ANDP",
    "MNChPD",
    "PCE",
    "TPR",
    "FPR"
  ),
  
  align = "c",
  
  booktabs = TRUE
) %>%
  
  kable_styling(
    
    full_width = FALSE,
    
    bootstrap_options = c(
      "striped",
      "hover",
      "condensed",
      "responsive"
    ),
    
    position = "center",
    
    font_size = 12
  ) %>%
  
  row_spec(
    0,
    bold = TRUE
  )

# =========================================================
# EXPORTAR
# =========================================================

write.csv(
  results_df,
  "metricas_completas.csv",
  row.names = FALSE
)

write.csv(
  summary_results,
  "resumen_metricas.csv",
  row.names = FALSE
)

# =========================================================
# TIEMPO TOTAL
# =========================================================

tiempo_fin <- Sys.time()

cat(
  "\nTiempo total de simulación:\n"
)

print(
  tiempo_fin - tiempo_inicio
)

