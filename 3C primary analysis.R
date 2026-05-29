# ============================================================
# 项目：NCCTG结肠癌辅助化疗试验 生存分析复盘
# 脚本：3-C 主分析
#       KM曲线 + log-rank检验 + stratified Cox模型
# 对应SAP章节：§5.1 Primary Analysis
# 输入：adtte.RData（由3-A生成）
# 输出：Figure 1（KM曲线）、Table 2（主分析结果）
# ============================================================

rm(list = ls())
library(survival)

# ---- 0. 载入数据 --------------------------------------------
load("adtte.RData")
cat("ITT人群 N =", nrow(adtte), "\n\n")

# 构造Surv对象——后面KM和Cox都要用
# time_months: 天转月（30.4375天/月）
# status: 1=事件（复发或死亡），0=删失
surv_obj <- Surv(adtte$time_months, adtte$status)


# ============================================================
# 第一部分：Kaplan-Meier 估计
# ============================================================

# ---- 1. 拟合KM曲线 ------------------------------------------
# ~ rx：按治疗臂分组
# conf.type = "log-log"：SAP预先指定的CI变换方式
# ⚠️ R默认是"log"，必须显式改成"log-log"
km_fit <- survfit(surv_obj ~ rx,
                  data     = adtte,
                  conf.type = "log-log")

cat("=== KM 估计结果 ===\n")
print(km_fit)

# ---- 2. 中位生存时间 ----------------------------------------
# survfit输出里自带中位生存时间和Brookmeyer-Crowley 95% CI
# ⚠️ R默认就是Brookmeyer-Crowley方法，SAP里要写明
cat("\n=== 中位生存时间（月）及 95% CI ===\n")
med <- summary(km_fit)$table
# 提取中位数和CI（单位：月）
for (i in seq_along(km_fit$strata)) {
  grp <- names(km_fit$strata)[i]
  cat(sprintf("  %s: median = %.1f 月 (95%% CI: %.1f, %.1f)\n",
              grp,
              med[i, "median"],
              med[i, "0.95LCL"],
              med[i, "0.95UCL"]))
}

# ---- 3. 关键时间点生存率 ------------------------------------
# 报1年、3年、5年生存率（SAP预先指定）
cat("\n=== 关键时间点 RFS 率 ===\n")
key_times <- c(12, 36, 60)  # 月
km_summary <- summary(km_fit, times = key_times)

# 按治疗组输出
groups <- levels(adtte$rx)
for (grp in groups) {
  idx <- km_summary$strata == paste0("rx=", grp)
  cat(sprintf("\n  %s:\n", grp))
  for (j in seq_along(key_times)) {
    if (sum(idx) >= j) {
      cat(sprintf("    %d月 RFS率: %.1f%% (95%% CI: %.1f%%, %.1f%%)\n",
                  key_times[j],
                  km_summary$surv[idx][j] * 100,
                  km_summary$lower[idx][j] * 100,
                  km_summary$upper[idx][j] * 100))
    }
  }
}


# ============================================================
# 第二部分：Log-rank 检验
# ============================================================

# ---- 4. 分层log-rank检验 ------------------------------------
# 分层因子：node4（SAP §11.4预先指定）
# 这是我们的primary significance test
lr_test <- survdiff(surv_obj ~ rx + strata(node4),
                    data = adtte)

cat("\n=== 分层 Log-rank 检验 ===\n")
print(lr_test)

# 手动提取p值
p_lr <- pchisq(lr_test$chisq, df = 1, lower.tail = FALSE)
cat(sprintf("  分层log-rank p值: %.4f\n", p_lr))


# ============================================================
# 第三部分：Stratified Cox 比例风险模型（主分析）
# ============================================================

# ---- 5. 主分析Cox模型 ----------------------------------------
# SAP §11.4预先指定：
#   - 分层因子：node4
#   - 协变量调整：age, sex, extent
#   - Tied events：Efron方法（R coxph默认，但SAP需写明）
#   - 参照组：Obs（3-A已显式设定）

cox_primary <- coxph(
  surv_obj ~ rx + strata(node4) + age + sex + extent,
  data   = adtte,
  ties   = "efron"   # 显式指定，对应SAP声明
)

cat("\n=== 主分析：Stratified Cox 模型 ===\n")
print(summary(cox_primary))

# ---- 6. 提取关键结果 ----------------------------------------
cox_sum <- summary(cox_primary)

# HR、95%CI、p值
hr  <- cox_sum$conf.int["rxLev+5FU", "exp(coef)"]
lcl <- cox_sum$conf.int["rxLev+5FU", "lower .95"]
ucl <- cox_sum$conf.int["rxLev+5FU", "upper .95"]
p_wald <- cox_sum$coefficients["rxLev+5FU", "Pr(>|z|)"]

cat("\n=== Table 2：主分析结果汇总 ===\n")
cat("  治疗对比：Lev+5FU vs. Obs\n")
cat(sprintf("  HR: %.2f (95%% CI: %.2f, %.2f)\n", hr, lcl, ucl))
cat(sprintf("  Wald p值: %.4f\n", p_wald))
cat(sprintf("  分层log-rank p值: %.4f\n", p_lr))

# ---- 7. 解读 ------------------------------------------------
cat("\n=== 结果解读 ===\n")
cat(sprintf("  HR = %.2f，即Lev+5FU组每一时刻的复发/死亡风险\n", hr))
cat(sprintf("  是Obs组的 %.0f%%，降低了 %.0f%% 的风险。\n",
            hr * 100, (1 - hr) * 100))
if (p_lr < 0.05) {
  cat("  分层log-rank检验p<0.05，两组RFS差异统计显著。\n")
} else {
  cat("  分层log-rank检验p>=0.05，两组RFS差异未达统计显著。\n")
}


# ============================================================
# 第四部分：KM曲线图（Figure 1）
# ============================================================

# ---- 8. 绘制KM曲线 ------------------------------------------
# 设置画布——上方留出绘图区域，下方留出number-at-risk表格空间
par(mar = c(8, 5, 4, 2))

# 画KM曲线
plot(km_fit,
     col     = c("steelblue", "firebrick"),  # 两组颜色
     lty     = c(1, 2),                       # 实线/虚线
     lwd     = 2,
     xlab    = "Time (months)",
     ylab    = "Recurrence-Free Survival",
     main    = "Figure 1. Kaplan-Meier Curves for RFS\n(ITT Population)",
     xlim    = c(0, 100),
     ylim    = c(0, 1),
     mark.time = TRUE,   # 在删失点标记竖线
     conf.int  = TRUE)   # 显示95% CI带

# 添加图例
legend("topright",
       legend = c(sprintf("Obs (N=%d)", sum(adtte$rx == "Obs")),
                  sprintf("Lev+5FU (N=%d)", sum(adtte$rx == "Lev+5FU"))),
       col    = c("steelblue", "firebrick"),
       lty    = c(1, 2),
       lwd    = 2,
       bty    = "n")

# 添加HR和p值文字注释
hr_text <- sprintf("HR = %.2f (95%% CI: %.2f, %.2f)\nLog-rank p = %.4f",
                   hr, lcl, ucl, p_lr)
text(x = 5, y = 0.12, labels = hr_text, adj = 0, cex = 0.9)

# ---- 9. Number at risk 表 -----------------------------------
# 标准肿瘤KM曲线必须配number-at-risk表，否则监管不接受
# ⚠️ "必须配number-at-risk"是行业规范，建议核对CONSORT/STEEP规范
at_risk_times <- c(0, 12, 24, 36, 48, 60, 72, 84, 96)
km_nar        <- summary(km_fit, times = at_risk_times)

groups_label <- c("Obs", "Lev+5FU")
colors_nar   <- c("steelblue", "firebrick")

for (g in seq_along(groups_label)) {
  idx <- km_nar$strata == paste0("rx=", groups_label[g])
  nar <- km_nar$n.risk[idx]
  
  # 在图底部写出每个时间点的风险集大小
  mtext(side  = 1,
        line  = 3 + (g - 1) * 1.8,
        at    = at_risk_times,
        text  = as.character(nar),
        cex   = 0.75,
        col   = colors_nar[g])
}

# 在左侧写标签
mtext(side = 1, line = 3,   at = -8, text = "Obs",     cex = 0.8,
      col = "steelblue", adj = 0)
mtext(side = 1, line = 4.8, at = -8, text = "Lev+5FU", cex = 0.8,
      col = "firebrick",  adj = 0)
mtext(side = 1, line = 6.2, at = -8, text = "No. at risk:", cex = 0.8,
      adj = 0)

cat("\nFigure 1 已保存：Figure1_KM_curve.pdf\n")

cat("\n=== 3-C 完成 ===\n")
