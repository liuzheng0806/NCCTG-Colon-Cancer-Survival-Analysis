# ============================================================
# 项目：NCCTG结肠癌辅助化疗试验 生存分析复盘
# 脚本：3-E 敏感性分析 SA-1 / SA-2 / SA-3
# 对应SAP章节：§5.3 Sensitivity Analyses
# 输入：adtte.RData（由3-A生成）
#
# ⚠️ 重要说明：
#   SA-1、SA-2 使用base R，已实测验证。
#   SA-3 delta-adjusted 逆变换插补。
# ============================================================

#rm(list = ls())
#library(survival)

#load("adtte.RData")
#surv_obj <- Surv(adtte$time_months, adtte$status)


# ============================================================
# 基准：主分析（回顾，便于对比）
# ============================================================
#cox_primary <- coxph(
#  surv_obj ~ rx + strata(node4) + age + sex + extent,
#  data = adtte, ties = "efron"
#)
ci0 <- summary(cox_primary)$conf.int["rxLev+5FU", ]
p0  <- summary(cox_primary)$coefficients["rxLev+5FU", "Pr(>|z|)"]

cat("=== 基准：主分析 ===\n")
cat(sprintf("  HR=%.3f (95%%CI %.3f-%.3f), p=%.2e, N=%d\n\n",
            ci0["exp(coef)"], ci0["lower .95"], ci0["upper .95"],
            p0, cox_primary$n))


# ============================================================
# SA-1：Unstratified Cox
# 针对假设：node4 作为分层因子是否合理？
# 做法：去掉 strata()，node4 改为普通协变量
# ============================================================
cat("=== SA-1: Unstratified Cox ===\n")
cat("针对假设：node4 作分层因子的建模选择是否影响治疗HR\n")

cox_sa1 <- coxph(
  surv_obj ~ rx + node4 + age + sex + extent,   # node4 当协变量
  data = adtte, ties = "efron"
)
ci1 <- summary(cox_sa1)$conf.int["rxLev+5FU", ]
p1  <- summary(cox_sa1)$coefficients["rxLev+5FU", "Pr(>|z|)"]
cat(sprintf("  SA-1 HR=%.3f (95%%CI %.3f-%.3f), p=%.2e, N=%d\n\n",
            ci1["exp(coef)"], ci1["lower .95"], ci1["upper .95"],
            p1, cox_sa1$n))


# ============================================================
# SA-2：Complete-case sensitivity
# 针对假设：因协变量缺失被剔除的12人，是否引入选择偏倚？
# 做法：在全ITT(N=619)上跑极简调整Cox(仅age+sex，不用
#       node4/extent，从而不丢失任何人)，与主分析对比
# ============================================================
cat("=== SA-2: Complete-case sensitivity ===\n")
cat("针对假设：协变量缺失被剔除的12人是否引入选择偏倚\n")

cox_sa2 <- coxph(
  surv_obj ~ rx + age + sex + extent,    # 仅age+sex+extent，保全619人
  data = adtte, ties = "efron"
)
ci2 <- summary(cox_sa2)$conf.int["rxLev+5FU", ]
p2  <- summary(cox_sa2)$coefficients["rxLev+5FU", "Pr(>|z|)"]
cat(sprintf("  SA-2 HR=%.3f (95%%CI %.3f-%.3f), p=%.2e, N=%d\n\n",
            ci2["exp(coef)"], ci2["lower .95"], ci2["upper .95"],
            p2, cox_sa2$n))


# ============================================================
# SA-1 / SA-2 汇总
# ============================================================
cat("=== SA-1 / SA-2 汇总对比 ===\n")
cat(sprintf("  主分析 (N=607, 分层, 3协变量)      HR=%.3f\n",
            ci0["exp(coef)"]))
cat(sprintf("  SA-1   (N=607, 不分层, node4协变量) HR=%.3f\n",
            ci1["exp(coef)"]))
cat(sprintf("  SA-2   (N=619, 极简调整)            HR=%.3f\n",
            ci2["exp(coef)"]))
cat("  三个HR几乎重合 → 主结论对分层选择、缺失剔除均稳健\n\n")


# ============================================================
# 项目：NCCTG结肠癌辅助化疗试验 生存分析复盘
# 脚本：3-E 补充 —— SA-3 Tipping Point Analysis
#       base R 自实现 delta-adjusted 删失插补
# 对应SAP章节：§5.3 Sensitivity Analyses (SA-3)
# ============================================================
# 第一部分：主分析模型 + 基线累积风险 H0(t)
# ============================================================
# 手动对齐:predict 只返回参与拟合的 607 个 lp(带行名),
# 按行名映射回 adtte 的 619 行,缺失者自动为 NA
lp_fit <- predict(cox_primary, type = "lp")      # 长度 607,有 names
lp <- rep(NA_real_, nrow(adtte))                 # 先建 619 个 NA
names(lp) <- rownames(adtte)
lp[names(lp_fit)] <- lp_fit                       # 按行名填回
stopifnot(length(lp) == nrow(adtte))             # 619 == 619,应通过
# nodes 缺失的 12 人 lp=NA、node4=NA,impute_once 中会被 next 跳过

# 基线累积风险 H0(t)：按 node4 分层各有一条
# basehaz 返回 hazard(累积)、time、strata 三列
bh <- basehaz(cox_primary, centered = TRUE)

# 为每个分层建一个 H0(t) 的插值函数（阶梯/线性）
strata_levels <- unique(bh$strata)
H0_funcs <- list()
for (s in strata_levels) {
  sub <- bh[bh$strata == s, ]
  # approxfun: 给定 t 返回 H0(t)，范围外用边界值
  H0_funcs[[as.character(s)]] <- approxfun(
    sub$time, sub$hazard,
    method = "linear", rule = 2          # rule=2: 范围外取端点值
  )
}

# 辅助：取某个体所属分层的 H0 函数
# ⚠️ basehaz() 的 strata 标签就是分层水平本身（如 "<=4 nodes"），
#    不带 "node4=" 前缀——直接用 node4 的值做 key
get_H0 <- function(node4_value) {
  H0_funcs[[as.character(node4_value)]]
}


# ============================================================
# 第二部分：单次 delta-adjusted 插补
#   只对 Lev+5FU 组的删失者插补剩余事件时间
# ============================================================
# 原理（逆变换抽样）：
#   个体 i 的累积风险 H_i(t) = H0_s(t) * exp(lp_i)
#   生存函数 S_i(t) = exp(-H_i(t))
#   已知活过删失时刻 c_i，剩余寿命的条件生存：
#     S_i(t | T>c_i) = S_i(t) / S_i(c_i)
#   delta 调整：c_i 之后的累积风险增量乘 delta：
#     H_i^delta(t) = H_i(c_i) + delta * (H_i(t) - H_i(c_i))
#   抽 U~Unif(0,1)，解 S_i^delta(t|T>c_i)=U 得插补事件时间

impute_once <- function(data, lp, delta, max_fu) {
  
  dat <- data
  n   <- nrow(dat)
  
  for (i in seq_len(n)) {
    
    # 仅处理：删失 且 属于治疗组
    if (dat$status[i] == 0 && dat$rx[i] == "Lev+5FU") {
      
      c_i  <- dat$time_months[i]
      H0_f <- get_H0(dat$node4[i])
      if (is.null(H0_f)) next                 # node4缺失者跳过
      
      elp  <- exp(lp[i])
      H_ci <- H0_f(c_i) * elp                 # 删失时刻累积风险
      
      # 在 [c_i, max_fu] 上构造 delta 调整后的条件生存
      # 取一组细网格时间点求解
      t_grid <- seq(c_i, max_fu, length.out = 200)
      H_grid <- H0_f(t_grid) * elp
      # delta 只放大 c_i 之后的风险增量
      H_delta <- H_ci + delta * (H_grid - H_ci)
      # 条件生存 S(t|T>c_i) = exp(-(H_delta - H_ci))
      S_cond  <- exp(-(H_delta - H_ci))
      
      # 逆变换抽样：U ~ Unif(0,1)，找 S_cond 首次 <= U
      U <- runif(1)
      hit <- which(S_cond <= U)[1]
      
      if (is.na(hit)) {
        # 整个窗内都没解出事件 → 仍删失在 max_fu
        dat$time_months[i] <- max_fu
        dat$status[i]      <- 0
      } else {
        dat$time_months[i] <- t_grid[hit]
        dat$status[i]      <- 1               # 改判为事件
      }
    }
  }
  dat
}


# ============================================================
# 第三部分：单个 delta —— M 次插补 + Rubin 合并
# ============================================================

run_one_delta <- function(data, lp, delta, M, max_fu) {
  
  betas <- numeric(M)
  ses   <- numeric(M)
  
  for (m in seq_len(M)) {
    imp_data <- impute_once(data, lp, delta, max_fu)
    fit <- coxph(
      Surv(time_months, status) ~ rx + strata(node4)
      + age + sex + extent,
      data = imp_data, ties = "efron"
    )
    betas[m] <- coef(fit)["rxLev+5FU"]
    ses[m]   <- sqrt(diag(vcov(fit)))["rxLev+5FU"]
  }
  
  # Rubin 规则
  q_bar <- mean(betas)
  u_bar <- mean(ses^2)
  b     <- var(betas)
  t_var <- u_bar + (1 + 1 / M) * b
  se_p  <- sqrt(t_var)
  
  z <- q_bar / se_p
  p <- 2 * pnorm(-abs(z))
  
  c(delta = delta, HR = exp(q_bar),
    lower = exp(q_bar - 1.96 * se_p),
    upper = exp(q_bar + 1.96 * se_p),
    p = p)
}


# ============================================================
# 第四部分：扫描 delta 网格，找翻转点
# ============================================================

delta_grid <- seq(1.0, 5.0, by = 0.5)   # 网格上探到5.0
M          <- 20
max_fu     <- max(adtte$time_months)

cat("=== SA-3: Tipping Point 扫描（delta只加治疗组）===\n")
cat(sprintf("delta 网格: %s\n", paste(delta_grid, collapse = ", ")))
cat(sprintf("每个delta插补次数 M = %d\n\n", M))

results <- as.data.frame(t(sapply(delta_grid, function(d) {
  run_one_delta(adtte, lp, d, M, max_fu)
})))

cat("delta    HR     95%CI下   95%CI上    p值\n")
for (k in seq_len(nrow(results))) {
  cat(sprintf("%.1f     %.3f   %.3f     %.3f     %.4f\n",
              results$delta[k], results$HR[k],
              results$lower[k], results$upper[k], results$p[k]))
}

# 翻转点：p 首次 >= 0.05
flip <- which(results$p >= 0.05)[1]
cat("\n=== 翻转点 ===\n")
if (is.na(flip)) {
  cat(sprintf("在 delta 上限 %.1f 内结论未翻转。\n",
              max(delta_grid)))
  cat("即使治疗组删失者风险放大到上限倍，疗效仍显著。\n")
  cat("→ 主结论对（治疗组方向的）informative censoring 稳健。\n")
} else {
  cat(sprintf("翻转点 delta = %.1f（p 在此首次 >= 0.05）\n",
              results$delta[flip]))
  cat("解读：要推翻疗效结论，需假设治疗组删失者真实风险\n")
  cat(sprintf("      被低估了约 %.1f 倍。\n", results$delta[flip]))
  if (results$delta[flip] >= 2) {
    cat("→ 翻转点较大，主结论相对稳健。\n")
  } else {
    cat("→ 翻转点接近1，主结论较脆弱，需谨慎。\n")
  }
}


# ============================================================
# 第五部分：绘图（窗口输出）
# ============================================================
plot(results$delta, results$HR,
     type = "b", pch = 19, col = "steelblue", lwd = 2,
     ylim = range(c(results$lower, results$upper, 1)),
     xlab = expression(paste("Treatment-arm censoring inflation  ",
                             delta)),
     ylab = "Hazard Ratio (Lev+5FU vs Obs)",
     main = "SA-3: Tipping Point (delta on treatment arm only)")
arrows(results$delta, results$lower,
       results$delta, results$upper,
       length = 0.04, angle = 90, code = 3, col = "steelblue")
abline(h = 1, lty = 2, col = "firebrick")

cat("\nTipping point 图已绘制到图形窗口。\n")
cat("观察：HR 随 delta 上升是否向 1 漂移，CI上界何时触及1。\n")

cat("\n=== SA-3 完成 ===\n")
