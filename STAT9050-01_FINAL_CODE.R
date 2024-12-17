# 1-(a) 단일 샘플에서 Cox 비례 위험 모형 적합 및 β 계수 추정

# 필요한 패키지 로드
library(MASS)       # 이변량 정규분포 생성을 위한 패키지
library(survival)   # 생존 분석(Cox 비례 위험 모형)을 위한 패키지

# 샘플 크기 및 파라미터 설정
set.seed(123)  # 재현성을 위한 시드 설정
n <- 30000     # 샘플 크기 설정
rho <- 0.5
gamma <- 1.5
beta <- c(log(2), log(2), -1)  # β1 = log(2), β2 = log(2), β3 = -1

# Z1, Z2, W1, W2 생성
# 상관계수 0.75를 가지는 이변량 정규분포로 Z1, Z2 생성
Sigma <- matrix(c(1, 0.75, 0.75, 1), nrow = 2)
Z <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
Z1 <- Z[, 1]
Z2 <- Z[, 2]

# W1은 표준 정규분포, W2는 베르누이 분포에서 생성
W1 <- rnorm(n)
W2 <- rbinom(n, 1, 0.5)

# 위험 함수에 따른 선형 예측값 및 생존 시간 T 생성
u <- runif(n)  # Uniform(0, 1) 분포에서 u 생성
lin_pred <- beta[1] * Z1 + beta[2] * W1 + beta[3] * W2  # 선형 예측값 계산

# Weibull 분포의 생존 시간 T 계산
# Weibull 생존 함수 S(t)에 따라 T 계산: T = ((-log(1 - u)) / (exp(lin_pred) * rho))^(1 / gamma)
T <- ((-log(1 - u)) / (exp(lin_pred) * rho))^(1 / gamma)

# Cox 비례 위험 모형 적합
# Surv() 함수로 생존 데이터를 정의하고 Cox 모델 적합
cox_model <- coxph(Surv(T) ~ Z1 + W1 + W2)  # Cox 모델 적합

# 모델 결과 요약 및 각 계수와 표준 오차 추정
summary_cox <- summary(cox_model)
beta_estimates <- summary_cox$coefficients[, "coef"]
se_estimates <- summary_cox$coefficients[, "se(coef)"]

# 결과 출력
cat(sprintf("추정된 β1: %.4f, 표준 오차: %.4f\n", beta_estimates[1], se_estimates[1]))
cat(sprintf("추정된 β2: %.4f, 표준 오차: %.4f\n", beta_estimates[2], se_estimates[2]))
cat(sprintf("추정된 β3: %.4f, 표준 오차: %.4f\n", beta_estimates[3], se_estimates[3]))


####################################

# 1-(b) 반복 시뮬레이션 및 요약 통계

# 반복 횟수 설정
iterations <- 100  # 반복 시뮬레이션 횟수
beta_results <- matrix(NA, nrow = iterations, ncol = 3)  # β 계수 추정치를 저장할 행렬
se_results <- matrix(NA, nrow = iterations, ncol = 3)    # 표준 오차 추정치를 저장할 행렬

for (i in seq_len(iterations)) {
  # 각 반복마다 새로운 Z1, Z2, W1, W2 및 T 생성
  
  # Z1과 Z2: 상관계수 0.75의 이변량 정규분포로부터 생성
  Z <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  
  # W1과 W2: 독립적인 표준 정규분포와 베르누이 분포로 생성
  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, 0.5)
  
  # Uniform(0,1)에서 샘플링하여 선형 예측값 및 생존 시간 T 생성
  u <- runif(n)
  lin_pred <- beta[1] * Z1 + beta[2] * W1 + beta[3] * W2
  T <- ((-log(1 - u)) / (exp(lin_pred) * rho))^(1 / gamma)
  
  # Cox 비례 위험 모형 적합
  cox_model <- coxph(Surv(T) ~ Z1 + W1 + W2)
  
  # 추정된 β 계수 및 표준 오차 저장
  beta_results[i, ] <- coef(cox_model)  # β 추정치 저장
  se_results[i, ] <- sqrt(diag(vcov(cox_model)))  # 표준 오차 저장
}

# 추정된 β 계수와 표준 오차의 평균 및 표준편차 계산
beta_means <- colMeans(beta_results, na.rm = TRUE)  # β 계수 평균
beta_sds <- apply(beta_results, 2, sd, na.rm = TRUE)  # β 계수 표준편차
se_means <- colMeans(se_results, na.rm = TRUE)  # 표준 오차 평균
se_sds <- apply(se_results, 2, sd, na.rm = TRUE)  # 표준 오차의 표준편차

# 결과 출력
cat("1-(b) 반복 시뮬레이션 결과:\n")
cat(sprintf("β1 평균: %.4f, 표준편차: %.4f\n", beta_means[1], beta_sds[1]))
cat(sprintf("β2 평균: %.4f, 표준편차: %.4f\n", beta_means[2], beta_sds[2]))
cat(sprintf("β3 평균: %.4f, 표준편차: %.4f\n", beta_means[3], beta_sds[3]))
cat(sprintf("β1 표준 오차 평균: %.4f, 표준편차: %.4f\n", se_means[1], se_sds[1]))
cat(sprintf("β2 표준 오차 평균: %.4f, 표준편차: %.4f\n", se_means[2], se_sds[2]))
cat(sprintf("β3 표준 오차 평균: %.4f, 표준편차: %.4f\n", se_means[3], se_sds[3]))


####################################

# 1-(c) 목표 검열율에 맞는 지수 분포의 λ 찾기 함수

# 필요한 패키지 로드
library(MASS)

# 목표 검열율에 맞는 λ 값을 찾기 위한 함수 정의
find_lambda <- function(target_rate, T) {
  # λ 값의 범위를 설정하고 목표 검열율에 가장 가까운 λ 값을 찾기
  result <- optimize(function(lambda) {
    C <- rexp(length(T), rate = lambda)  # 검열 시간 C 생성
    censor_rate <- mean(C < T)           # 실제 검열율 계산
    return(abs(censor_rate - target_rate)) # 목표 검열율과 차이를 최소화
  }, interval = c(1e-5, 10))              # λ 검색 구간 설정
  
  return(result$minimum)                  # 최적화된 λ 반환
}

# 샘플 크기 및 파라미터 설정
set.seed(123)
large_n <- 1e6
rho <- 0.5
gamma <- 1.5
beta <- c(log(2), log(2), -1)

# 큰 샘플 크기에서 Z1, Z2, W1, W2 및 생존 시간 T 생성
Sigma <- matrix(c(1, 0.75, 0.75, 1), nrow = 2)
Z <- mvrnorm(large_n, mu = c(0, 0), Sigma = Sigma)
W1 <- rnorm(large_n)
W2 <- rbinom(large_n, 1, 0.5)
u <- runif(large_n)
lin_pred <- beta[1] * Z[,1] + beta[2] * W1 + beta[3] * W2
T <- ((-log(1 - u)) / (exp(lin_pred) * rho))^(1 / gamma)

# 목표 검열율 설정
censor_rates <- c(0.10, 0.30, 0.90, 0.95, 0.99)

# 각 검열율에 맞는 λ 값 찾기
lambda_values <- sapply(censor_rates, find_lambda, T = T)

# 결과 출력
cat("1-(c) 목표 검열율에 따른 λ 값:\n")
for (i in seq_along(censor_rates)) {
  cat(sprintf("검열율 %.2f%%에 대한 λ: %.4f\n", censor_rates[i] * 100, lambda_values[i]))
}
####################################

# 1-(d) 각 검열율에 따른 반복 시뮬레이션

results <- lapply(seq_along(censor_rates), function(rate_idx) {
  lambda <- lambda_values[rate_idx]  # 각 검열율에 맞는 λ 값 가져오기
  
  # λ가 NA일 경우 해당 검열율에 대한 시뮬레이션 건너뜀
  if (is.na(lambda)) {
    warning(sprintf("검열율 %.2f%% 에 대한 λ 값이 없습니다. 시뮬레이션을 건너뜁니다.", censor_rates[rate_idx] * 100))
    return(NULL)
  }
  
  beta_matrix <- matrix(NA, nrow = iterations, ncol = 3)  # β 계수 저장 행렬
  se_matrix <- matrix(NA, nrow = iterations, ncol = 3)    # 표준 오차 저장 행렬
  
  for (i in seq_len(iterations)) {
    # Z1, Z2, W1, W2 및 T 생성
    Z <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
    Z1 <- Z[, 1]
    Z2 <- Z[, 2]
    W1 <- rnorm(n)
    W2 <- rbinom(n, 1, 0.5)
    u <- runif(n)
    lin_pred <- beta[1] * Z1 + beta[2] * W1 + beta[3] * W2
    T <- ((-log(1 - u)) / (exp(lin_pred) * rho))^(1 / gamma)
    
    # 검열 시간 C 생성
    C <- rexp(n, rate = lambda)
    
    # 관측 시간 X와 이벤트 인디케이터 Delta 생성
    X <- pmin(T, C)  # T와 C의 최소값을 관측 시간 X로 설정
    Delta <- as.numeric(T <= C)  # 이벤트 여부를 Delta로 설정 (1: 이벤트 발생, 0: 검열)
    
    # 유효한 관측값이 있는지 확인 (모든 관측값이 검열되면 모델 적합 불가)
    if (all(Delta == 0)) {
      warning("모든 관측값이 검열되었습니다. 모델 적합이 불가능합니다.")
      next  # 다음 반복으로 건너뜀
    }
    
    # Cox 비례 위험 모형 적합
    cox_model <- coxph(Surv(X, Delta) ~ Z1 + W1 + W2)
    
    # 추정된 β 계수 및 표준 오차 저장
    beta_matrix[i, ] <- coef(cox_model)
    se_matrix[i, ] <- sqrt(diag(vcov(cox_model)))
  }
  
  list(
    beta_means = colMeans(beta_matrix, na.rm = TRUE),
    beta_sds = apply(beta_matrix, 2, sd, na.rm = TRUE),
    se_means = colMeans(se_matrix, na.rm = TRUE),
    se_sds = apply(se_matrix, 2, sd, na.rm = TRUE)
  )
})

# 각 검열율에 대한 결과 출력
cat("1-(d) 검열율별 반복 시뮬레이션 결과:\n")
for (i in seq_along(results)) {
  res <- results[[i]]
  
  # 시뮬레이션이 NULL인 경우 건너뜀
  if (is.null(res)) next
  
  cat(sprintf("검열율 %.2f%% 결과:\n", censor_rates[i] * 100))
  cat(sprintf("β1 평균: %.4f, 표준편차: %.4f\n", res$beta_means[1], res$beta_sds[1]))
  cat(sprintf("β2 평균: %.4f, 표준편차: %.4f\n", res$beta_means[2], res$beta_sds[2]))
  cat(sprintf("β3 평균: %.4f, 표준편차: %.4f\n", res$beta_means[3], res$beta_sds[3]))
  cat(sprintf("β1 표준 오차 평균: %.4f, 표준편차: %.4f\n", res$se_means[1], res$se_sds[1]))
  cat(sprintf("β2 표준 오차 평균: %.4f, 표준편차: %.4f\n", res$se_means[2], res$se_sds[2]))
  cat(sprintf("β3 표준 오차 평균: %.4f, 표준편차: %.4f\n", res$se_means[3], res$se_sds[3]))
}


####################################
####################################


# 2-(a) CASE COHORT
library(MASS)
library(survival)

# 설정 값
n <- 30000
iterations <- 500
beta1 <- log(2)
beta2 <- log(2)
beta3 <- -1
rho <- 0.5
gamma <- 1.5
censoring_rate <- 0.99

set.seed(123)

# 데이터 생성 함수
generate_data <- function() {
  Sigma <- matrix(c(1, 0.75, 0.75, 1), 2, 2)
  Z <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, 0.5)
  
  # 생존 시간 생성
  linear_pred <- beta1 * Z1 + beta2 * W1 + beta3 * W2
  U <- runif(n)
  T <- (-log(1 - U) / (rho * exp(linear_pred)))^(1 / gamma)
  
  # 검열 시간 생성
  lambda_censor <- quantile(T, probs = 1 - censoring_rate)
  C <- rexp(n, rate = 1 / lambda_censor)
  
  # 관측 시간 및 상태
  X <- pmin(T, C)
  delta <- as.numeric(T <= C)
  
  data.frame(X = X, delta = delta, Z1 = Z1, Z2 = Z2, W1 = W1, W2 = W2)
}

# Subcohort sizes
subcohort_sizes <- c(100, 300)
all_results <- list()

# Loop over subcohort sizes
for (subcohort_size in subcohort_sizes) {
  results <- vector("list", iterations)  # 결과 저장 리스트 초기화
  
  for (iter in 1:iterations) {
    data <- generate_data()  # 각 반복에서 데이터 생성
    
    # 실패 사례와 대조군 분리
    failures <- which(data$delta == 1)
    non_failures <- which(data$delta == 0)
    
    # 무작위 서브코호트 구성
    random_subcohort <- sample(1:n, subcohort_size, replace = FALSE)
    
    # 실패 사례 추가
    unique_failures <- setdiff(failures, random_subcohort)  # 중복되지 않은 실패 사례
    final_subcohort_indices <- c(random_subcohort, unique_failures)
    
    data_sub <- data[final_subcohort_indices, ]
    
    # 가중치 계산
    n_censored <- sum(data_sub$delta == 0)
    if (n_censored > 0) {
      weights <- ifelse(data_sub$delta == 1, 1, n / n_censored)
    } else {
      next  # 대조군이 없으면 다음 반복으로 이동
    }
    
    # Cox 모델 적합
    fit <- tryCatch({
      coxph(Surv(X, delta) ~ Z1 + W1 + W2, data = data_sub, weights = weights)
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      beta <- coef(fit)
      se <- sqrt(diag(vcov(fit)))
      results[[iter]] <- list(beta = beta, se = se)
    }
  }
  
  # 결과 요약
  valid_results <- results[!sapply(results, is.null)]  # 유효한 결과만 추출
  
  if (length(valid_results) > 0) {
    beta_estimates <- do.call(rbind, lapply(valid_results, `[[`, "beta"))
    se_estimates <- do.call(rbind, lapply(valid_results, `[[`, "se"))
    
    # 평균 및 표준편차 계산
    beta_means <- colMeans(beta_estimates, na.rm = TRUE)
    beta_sds <- apply(beta_estimates, 2, sd, na.rm = TRUE)
    se_means <- colMeans(se_estimates, na.rm = TRUE)
    se_sds <- apply(se_estimates, 2, sd, na.rm = TRUE)
    
    all_results[[as.character(subcohort_size)]] <- list(
      beta_means = beta_means,
      beta_sds = beta_sds,
      se_means = se_means,
      se_sds = se_sds
    )
  }
}

# 결과 출력
for (subcohort_size in subcohort_sizes) {
  cat(sprintf("\nCase-Cohort Results (Subcohort Size = %d):\n", subcohort_size))
  if (!is.null(all_results[[as.character(subcohort_size)]])) {
    cat("Beta 추정치 평균:", all_results[[as.character(subcohort_size)]]$beta_means, "\n")
    cat("Beta 추정치 표준편차:", all_results[[as.character(subcohort_size)]]$beta_sds, "\n")
    cat("표준 오차 평균:", all_results[[as.character(subcohort_size)]]$se_means, "\n")
    cat("표준 오차 표준편차:", all_results[[as.character(subcohort_size)]]$se_sds, "\n")
  } else {
    cat("결과 없음\n")
  }
}

####################################################

# 2-(b) Nested Case-Control 

library(MASS)
library(survival)

# 설정 값
n <- 30000
control_sizes <- c(1, 5)  # 대조군 수 설정: 1명, 5명
iterations <- 3000
beta1 <- log(2)
beta2 <- log(2)
beta3 <- -1
rho <- 0.5
gamma <- 1.5
censoring_rate <- 0.99

set.seed(123)

# 데이터 생성 함수
generate_data <- function() {
  Sigma <- matrix(c(1, 0.75, 0.75, 1), 2, 2)
  Z <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, 0.5)
  
  # 생존 시간 생성
  linear_pred <- beta1 * Z1 + beta2 * W1 + beta3 * W2
  U <- runif(n)
  T <- (-log(1 - U) / (rho * exp(linear_pred)))^(1 / gamma)
  
  # 검열 시간 생성
  lambda_censor <- quantile(T, probs = 1 - censoring_rate)
  C <- rexp(n, rate = 1 / lambda_censor)
  
  # 관측 시간 및 상태
  X <- pmin(T, C)
  delta <- as.numeric(T <= C)
  
  data.frame(X = X, delta = delta, Z1 = Z1, Z2 = Z2, W1 = W1, W2 = W2)
}

# Nested Case-Control 분석
nested_results <- list()

for (control_size in control_sizes) {
  results <- vector("list", iterations)  # 결과 저장 리스트 초기화
  
  for (iter in 1:iterations) {
    data <- generate_data()  # 각 반복에서 데이터 생성
    
    # 사건 발생자(케이스) 확인
    cases <- which(data$delta == 1)
    
    sampled_indices <- numeric()  # Nested Case-Control 서브셋 인덱스 저장
    
    for (case_id in cases) {
      event_time <- data$X[case_id]
      
      # 위험 집단: 사건 발생 시점까지 생존한 피험자
      risk_set <- which(data$X > event_time)
      
      if (length(risk_set) < control_size) {
        # 위험 집단 크기가 컨트롤 수보다 작으면 전체를 포함
        sampled_case_control <- c(case_id, risk_set)
      } else {
        # 위험 집단에서 컨트롤 샘플링
        controls <- sample(risk_set, control_size, replace = FALSE)
        sampled_case_control <- c(case_id, controls)
      }
      
      sampled_indices <- c(sampled_indices, sampled_case_control)
    }
    
    # Nested Case-Control 데이터 구성
    data_sub <- data[unique(sampled_indices), ]
    
    # Cox 모델 적합
    fit <- tryCatch({
      coxph(Surv(X, delta) ~ Z1 + W1 + W2, data = data_sub)
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      beta <- coef(fit)
      se <- sqrt(diag(vcov(fit)))
      results[[iter]] <- list(beta = beta, se = se)
    }
  }
  
  # 결과 요약
  valid_results <- results[!sapply(results, is.null)]  # 유효한 결과만 추출
  
  if (length(valid_results) > 0) {
    beta_estimates <- do.call(rbind, lapply(valid_results, `[[`, "beta"))
    se_estimates <- do.call(rbind, lapply(valid_results, `[[`, "se"))
    
    # 평균 및 표준편차 계산
    beta_means <- colMeans(beta_estimates, na.rm = TRUE)
    beta_sds <- apply(beta_estimates, 2, sd, na.rm = TRUE)
    se_means <- colMeans(se_estimates, na.rm = TRUE)
    se_sds <- apply(se_estimates, 2, sd, na.rm = TRUE)
    
    # 결과 저장
    nested_results[[as.character(control_size)]] <- list(
      beta_means = beta_means,
      beta_sds = beta_sds,
      se_means = se_means,
      se_sds = se_sds
    )
  }
}

# 결과 출력
for (control_size in control_sizes) {
  cat(sprintf("Nested Case-Control (Control Size = %d):\n", control_size))
  if (!is.null(nested_results[[as.character(control_size)]])) {
    cat("Beta 추정치 평균:", nested_results[[as.character(control_size)]]$beta_means, "\n")
    cat("Beta 추정치 표준편차:", nested_results[[as.character(control_size)]]$beta_sds, "\n")
    cat("표준 오차 평균:", nested_results[[as.character(control_size)]]$se_means, "\n")
    cat("표준 오차 표준편차:", nested_results[[as.character(control_size)]]$se_sds, "\n\n")
  } else {
    cat("결과 없음\n")
  }
}

####################################
####################################

# 3-(a),(b),(c)
# 필요한 패키지 로드
library(MASS)
library(survival)
library(dplyr)

# 설정 값
n <- 30000
beta1 <- log(2)
beta2 <- log(2)
beta3 <- -1
rho <- 0.5
gamma <- 1.5
censoring_rate <- 0.99

set.seed(123)

# 데이터 생성 함수
generate_data <- function() {
  Sigma <- matrix(c(1, 0.75, 0.75, 1), 2, 2)
  Z <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, 0.5)
  
  linear_pred <- beta1 * Z1 + beta2 * W1 + beta3 * W2
  U <- runif(n)
  T <- (-log(1 - U) / (rho * exp(linear_pred)))^(1 / gamma)
  
  lambda_censor <- quantile(T, probs = 1 - censoring_rate)
  C <- rexp(n, rate = 1 / lambda_censor)
  
  X <- pmin(T, C)
  delta <- as.numeric(T <= C)
  
  data.frame(X = X, delta = delta, Z1 = Z1, Z2 = Z2, W1 = W1, W2 = W2)
}

# 서브코호트 구성 함수
create_subcohort <- function(data, subcohort_size, stratified = FALSE) {
  failures <- which(data$delta == 1)
  non_failures <- which(data$delta == 0)
  
  if (stratified) {
    fit <- coxph(Surv(X, delta) ~ W1 + W2, data = data)
    influence_scores <- residuals(fit, type = "dfbeta") %>% rowSums() %>% `^`(2)
    influence_scores <- pmax(influence_scores, 1e-5)
    influence_weights <- exp(influence_scores) / sum(exp(influence_scores))
    
    data$quantile_group <- ntile(data$Z2, 5)
    stratified_sample <- unlist(lapply(split(non_failures, data$quantile_group[non_failures]), function(group) {
      sample(group, size = round(subcohort_size / 5), prob = influence_weights[group], replace = FALSE)
    }))
    
    final_subcohort_indices <- c(stratified_sample, failures)
  } else {
    random_subcohort <- sample(1:n, subcohort_size, replace = FALSE)
    unique_failures <- setdiff(failures, random_subcohort)
    final_subcohort_indices <- c(random_subcohort, unique_failures)
  }
  
  data[final_subcohort_indices, ]
}

# 분석 함수
analyze_subcohort <- function(data_sub, n) {
  n_censored <- sum(data_sub$delta == 0)
  if (n_censored > 0) {
    weights <- ifelse(data_sub$delta == 1, 1, n / n_censored)
  } else {
    return(NULL)
  }
  
  fit <- tryCatch({
    coxph(Surv(X, delta) ~ Z1 + W1 + W2, data = data_sub, weights = weights)
  }, error = function(e) NULL)
  
  if (!is.null(fit)) {
    beta <- coef(fit)
    se <- sqrt(diag(vcov(fit)))
    return(list(beta = beta, se = se))
  }
  return(NULL)
}

# 실행 함수
run_case_cohort_analysis <- function(subcohort_size, iterations, stratified = FALSE) {
  results <- vector("list", iterations)
  
  for (iter in 1:iterations) {
    data <- generate_data()
    data_sub <- create_subcohort(data, subcohort_size, stratified)
    result <- analyze_subcohort(data_sub, n)
    if (!is.null(result)) {
      results[[iter]] <- result
    }
  }
  
  valid_results <- results[!sapply(results, is.null)]
  if (length(valid_results) > 0) {
    beta_estimates <- do.call(rbind, lapply(valid_results, `[[`, "beta"))
    se_estimates <- do.call(rbind, lapply(valid_results, `[[`, "se"))
    
    beta_means <- colMeans(beta_estimates, na.rm = TRUE)
    beta_sds <- apply(beta_estimates, 2, sd, na.rm = TRUE)
    se_means <- colMeans(se_estimates, na.rm = TRUE)
    se_sds <- apply(se_estimates, 2, sd, na.rm = TRUE)
    
    return(list(beta_means = beta_means, beta_sds = beta_sds, se_means = se_means, se_sds = se_sds))
  }
  return(NULL)
}

# 통합 비교 함수
compare_methods <- function(subcohort_size, iterations) {
  original_results <- run_case_cohort_analysis(subcohort_size, iterations, stratified = FALSE)
  calibrated_results <- run_case_cohort_analysis(subcohort_size, iterations, stratified = TRUE)
  
  if (!is.null(original_results) & !is.null(calibrated_results)) {
    comparison <- tibble(
      Method = c("Original", "Calibrated"),
      Beta1_Bias = c(original_results$beta_means[1] - log(2), calibrated_results$beta_means[1] - log(2)),
      Beta2_Bias = c(original_results$beta_means[2] - log(2), calibrated_results$beta_means[2] - log(2)),
      Beta3_Bias = c(original_results$beta_means[3] - (-1), calibrated_results$beta_means[3] - (-1)),
      Beta1_SE = c(original_results$se_means[1], calibrated_results$se_means[1]),
      Beta2_SE = c(original_results$se_means[2], calibrated_results$se_means[2]),
      Beta3_SE = c(original_results$se_means[3], calibrated_results$se_means[3])
    )
    return(comparison)
  } else {
    return(NULL)
  }
}

# 실행 및 결과 비교
subcohort_size <- 200
iterations <- 100
comparison <- compare_methods(subcohort_size, iterations)
print(comparison)

#####################################

# iteration, subcohort size 다양하게 simulation

library(MASS)
library(survival)
library(dplyr)

# 설정 값
n <- 30000
beta1 <- log(2)
beta2 <- log(2)
beta3 <- -1
rho <- 0.5
gamma <- 1.5
censoring_rate <- 0.99

set.seed(123)

# 데이터 생성 함수
generate_data <- function() {
  Sigma <- matrix(c(1, 0.75, 0.75, 1), 2, 2)
  Z <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, 0.5)
  
  linear_pred <- beta1 * Z1 + beta2 * W1 + beta3 * W2
  U <- runif(n)
  T <- (-log(1 - U) / (rho * exp(linear_pred)))^(1 / gamma)
  
  lambda_censor <- quantile(T, probs = 1 - censoring_rate)
  C <- rexp(n, rate = 1 / lambda_censor)
  
  X <- pmin(T, C)
  delta <- as.numeric(T <= C)
  
  data.frame(X = X, delta = delta, Z1 = Z1, Z2 = Z2, W1 = W1, W2 = W2)
}

# 서브코호트 구성 함수
create_subcohort <- function(data, subcohort_size, stratified = FALSE) {
  failures <- which(data$delta == 1)
  non_failures <- which(data$delta == 0)
  
  if (stratified) {
    fit <- coxph(Surv(X, delta) ~ W1 + W2, data = data)
    influence_scores <- residuals(fit, type = "dfbeta") %>% rowSums() %>% `^`(2)
    influence_scores <- pmax(influence_scores, 1e-5)
    influence_weights <- exp(influence_scores) / sum(exp(influence_scores))
    
    data$quantile_group <- ntile(data$Z2, 5)
    stratified_sample <- unlist(lapply(split(non_failures, data$quantile_group[non_failures]), function(group) {
      sample(group, size = round(subcohort_size / 5), prob = influence_weights[group], replace = FALSE)
    }))
    
    final_subcohort_indices <- c(stratified_sample, failures)
  } else {
    random_subcohort <- sample(1:n, subcohort_size, replace = FALSE)
    unique_failures <- setdiff(failures, random_subcohort)
    final_subcohort_indices <- c(random_subcohort, unique_failures)
  }
  
  data[final_subcohort_indices, ]
}

# 분석 함수
analyze_subcohort <- function(data_sub, n) {
  n_censored <- sum(data_sub$delta == 0)
  if (n_censored > 0) {
    weights <- ifelse(data_sub$delta == 1, 1, n / n_censored)
  } else {
    return(NULL)
  }
  
  fit <- tryCatch({
    coxph(Surv(X, delta) ~ Z1 + W1 + W2, data = data_sub, weights = weights)
  }, error = function(e) NULL)
  
  if (!is.null(fit)) {
    beta <- coef(fit)
    se <- sqrt(diag(vcov(fit)))
    return(list(beta = beta, se = se))
  }
  return(NULL)
}

# 실행 함수
run_case_cohort_analysis <- function(subcohort_size, iterations, stratified = FALSE) {
  results <- vector("list", iterations)
  
  for (iter in 1:iterations) {
    data <- generate_data()
    data_sub <- create_subcohort(data, subcohort_size, stratified)
    result <- analyze_subcohort(data_sub, n)
    if (!is.null(result)) {
      results[[iter]] <- result
    }
  }
  
  valid_results <- results[!sapply(results, is.null)]
  if (length(valid_results) > 0) {
    beta_estimates <- do.call(rbind, lapply(valid_results, `[[`, "beta"))
    se_estimates <- do.call(rbind, lapply(valid_results, `[[`, "se"))
    
    beta_means <- colMeans(beta_estimates, na.rm = TRUE)
    beta_sds <- apply(beta_estimates, 2, sd, na.rm = TRUE)
    se_means <- colMeans(se_estimates, na.rm = TRUE)
    se_sds <- apply(se_estimates, 2, sd, na.rm = TRUE)
    
    return(list(beta_means = beta_means, beta_sds = beta_sds, se_means = se_means, se_sds = se_sds))
  }
  return(NULL)
}

# 통합 비교 함수
compare_methods <- function(subcohort_size, iterations) {
  original_results <- run_case_cohort_analysis(subcohort_size, iterations, stratified = FALSE)
  calibrated_results <- run_case_cohort_analysis(subcohort_size, iterations, stratified = TRUE)
  
  if (!is.null(original_results) & !is.null(calibrated_results)) {
    comparison <- tibble(
      Subcohort_Size = subcohort_size,
      Iterations = iterations,
      Method = c("Original", "Calibrated"),
      Beta1_Bias = c(original_results$beta_means[1] - log(2), calibrated_results$beta_means[1] - log(2)),
      Beta2_Bias = c(original_results$beta_means[2] - log(2), calibrated_results$beta_means[2] - log(2)),
      Beta3_Bias = c(original_results$beta_means[3] - (-1), calibrated_results$beta_means[3] - (-1)),
      Beta1_SE = c(original_results$se_means[1], calibrated_results$se_means[1]),
      Beta2_SE = c(original_results$se_means[2], calibrated_results$se_means[2]),
      Beta3_SE = c(original_results$se_means[3], calibrated_results$se_means[3])
    )
    return(comparison)
  } else {
    return(NULL)
  }
}

# 모든 케이스 실행 및 결과 통합
run_all_cases <- function() {
  subcohort_sizes <- c(300, 500, 1000)
  iteration_counts <- c(1000, 3000)
  
  results <- lapply(subcohort_sizes, function(size) {
    lapply(iteration_counts, function(iters) {
      compare_methods(size, iters)
    })
  })
  
  # 결과를 데이터프레임으로 통합
  do.call(rbind, lapply(results, function(res) do.call(rbind, res)))
}

# 결과 실행
final_results <- run_all_cases()
print(final_results)


####################################

# 3-(d)
library(survival)
library(dplyr)
library(tidyr)
library(kableExtra)

# Real Data
data <- read.csv("C:/Users/gimyo/Desktop/mimic3_final.csv")
colnames(data) <- c("T", "Z1", "HeartRate", "Height", "Z2", 
                    "OxygenSaturation", "RespRate", "Temperature", 
                    "Weight", "EyeOpening", "VerbalResponse", "delta")

# Original 방법 실행 함수
run_original <- function(data, subcohort_size) {
  set.seed(123)
  
  failures <- which(data$delta == 1)
  non_failures <- which(data$delta == 0)
  random_subcohort <- sample(non_failures, subcohort_size, replace = FALSE)
  unique_failures <- setdiff(failures, random_subcohort)
  final_subcohort_indices <- c(random_subcohort, unique_failures)
  data_sub <- data[final_subcohort_indices, ]
  
  n_censored <- sum(data_sub$delta == 0)
  weights <- ifelse(data_sub$delta == 1, 1, nrow(data) / n_censored)
  
  tryCatch({
    fit <- coxph(Surv(T, delta) ~ Z1 + HeartRate + Height + OxygenSaturation + RespRate + Temperature + Weight + EyeOpening + VerbalResponse, 
                 data = data_sub, weights = weights)
    list(beta = coef(fit), se = sqrt(diag(vcov(fit))))
  }, error = function(e) NULL)
}

# Calibrated 방법 실행 함수
run_calibrated <- function(data, subcohort_size) {
  set.seed(123)
  
  failures <- which(data$delta == 1)
  non_failures <- which(data$delta == 0)
  data$quantile_group <- ntile(data$Z2, 5)
  fit <- coxph(Surv(T, delta) ~ HeartRate + Height, data = data)
  influence_scores <- residuals(fit, type = "dfbeta") %>% rowSums() %>% `^`(2)
  influence_weights <- pmax(influence_scores, 1e-5)
  influence_weights <- exp(influence_weights) / sum(exp(influence_weights))
  
  stratified_sample <- unlist(lapply(split(non_failures, data$quantile_group[non_failures]), function(group) {
    sample(group, size = round(subcohort_size / 5), prob = influence_weights[group], replace = FALSE)
  }))
  
  final_subcohort_indices <- c(stratified_sample, failures)
  data_sub <- data[final_subcohort_indices, ]
  
  n_censored <- sum(data_sub$delta == 0)
  weights <- ifelse(data_sub$delta == 1, 1, nrow(data) / n_censored)
  
  tryCatch({
    fit <- coxph(Surv(T, delta) ~ Z1 + HeartRate + Height + OxygenSaturation + RespRate + Temperature + Weight + EyeOpening + VerbalResponse, 
                 data = data_sub, weights = weights)
    list(beta = coef(fit), se = sqrt(diag(vcov(fit))))
  }, error = function(e) NULL)
}

# Gold Standard 실행 함수
run_gold_standard <- function(data) {
  tryCatch({
    fit <- coxph(Surv(T, delta) ~ Z1 + HeartRate + Height + OxygenSaturation + RespRate + Temperature + Weight + EyeOpening + VerbalResponse, data = data)
    list(beta = coef(fit), se = sqrt(diag(vcov(fit))))
  }, error = function(e) NULL)
}

# 통합 비교 함수
compare_methods_all_variables <- function(data, subcohort_size) {
  original <- run_original(data, subcohort_size)
  calibrated <- run_calibrated(data, subcohort_size)
  gold_standard <- run_gold_standard(data)
  
  if (!is.null(original) & !is.null(calibrated) & !is.null(gold_standard)) {
    variables <- c("Z1", "HeartRate", "Height", "OxygenSaturation", 
                   "RespRate", "Temperature", "Weight", 
                   "EyeOpening", "VerbalResponse")
    
    tibble(
      Subcohort_Size = subcohort_size,
      Variable = variables,
      `Original Beta(SE)` = sapply(variables, function(var) {
        paste0(round(original$beta[var], 4), " (", round(original$se[var], 4), ")")
      }),
      `Calibrated Beta(SE)` = sapply(variables, function(var) {
        paste0(round(calibrated$beta[var], 4), " (", round(calibrated$se[var], 4), ")")
      }),
      `Gold Standard Beta(SE)` = sapply(variables, function(var) {
        paste0(round(gold_standard$beta[var], 4), " (", round(gold_standard$se[var], 4), ")")
      })
    )
  } else {
    NULL
  }
}

# 모든 서브코호트 크기에 대해 실행
run_all_subcohort_sizes <- function(data, subcohort_sizes) {
  results <- lapply(subcohort_sizes, function(size) {
    compare_methods_all_variables(data, size)
  })
  
  do.call(rbind, results)
}

# 실행 및 결과
subcohort_sizes <- c(1000, 2000, 3000)
comparison_results_all <- run_all_subcohort_sizes(data, subcohort_sizes)

# 결과 출력
if (!is.null(comparison_results_all)) {
  comparison_results_all %>%
    kable() %>%
    kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
} else {
  print("One or more methods failed.")
}

####################################
############### END ################
####################################
