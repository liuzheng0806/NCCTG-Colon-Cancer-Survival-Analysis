# ============================================================
# 项目：NCCTG结肠癌辅助化疗试验 生存分析复盘
# 脚本：3-A 环境搭建 + 数据加载 + ITT人群构造
# 数据来源：survival::colon (R package)
# 对应SAP章节：§3 Study Design, §4 Analysis Populations
# ============================================================

# ---- 0. 环境准备 --------------------------------------------

# 清空工作环境，保证可重复性
rm(list = ls())

# 记录R版本（对应SAP §9 Software & Reproducibility）
# 在真实SAP里这个版本号要写死
cat("R version:\n")
print(R.version.string)

# 加载必要包
# 只依赖 survival（base R），避免额外依赖
library(survival)

cat("\nsurvival 包版本:", as.character(packageVersion("survival")), "\n")


# ---- 1. 加载原始数据 ----------------------------------------

# colon数据集说明：
#   - 每个患者有2行，由etype区分：
#     etype = 1: recurrence（复发事件）→ 我们用这个构造RFS
#     etype = 2: death（死亡事件）      → 本项目不用，留作OS备用
#   - 治疗组rx：Obs / Lev / Lev+5FU（三臂）

data(cancer, package = "survival")  # colon在cancer里

# 核查原始数据结构
cat("\n=== 原始数据结构 ===\n")
cat("总行数（应为929患者×2行 = 1858）:", nrow(colon), "\n")
cat("变量列表:\n")
print(names(colon))

cat("\n治疗组×etype 分布（原始三臂）:\n")
print(table(colon$rx, colon$etype,
            dnn = c("治疗组", "etype(1=recurrence,2=death)")))


# ---- 2. 构造ITT分析数据集 -----------------------------------
# 对应SAP §4：ITT Population
# 规则：
#   (a) 只保留recurrence事件行（etype == 1）
#   (b) 只保留Obs和Lev+5FU两个臂（剔除Lev单药组）
#   (c) 按SAP预先指定，Obs为参照组

# (a) + (b) 过滤
adtte <- colon[colon$etype == 1 &
                 colon$rx %in% c("Obs", "Lev+5FU"), ]

# (c) 显式设定因子水平和参照组
#     第一个水平 = Cox参照组
#     Obs在前 → HR解读为"Lev+5FU相对于Obs的风险比"
adtte$rx <- factor(adtte$rx, levels = c("Obs", "Lev+5FU"))

# 派生变量1：node4（是否>4个阳性淋巴结）
# 用途：SAP预设的分层因子
adtte$node4 <- factor(
  ifelse(adtte$nodes > 4, ">4 nodes", "<=4 nodes"),
  levels = c("<=4 nodes", ">4 nodes")
)

# 派生变量2：时间单位转换 天→月
# 用途：报中位生存时间时用月更直观（临床习惯）
# 转换系数：365.25/12 ≈ 30.4375天/月
adtte$time_months <- adtte$time / (365.25 / 12)

# 派生变量3：differ转因子（肿瘤分化程度）
# 原始编码：1=well, 2=moderate, 3=poor
adtte$differ <- factor(adtte$differ,
                       levels = c(1, 2, 3),
                       labels = c("Well", "Moderate", "Poor"))

# 派生变量4：extent转因子（肿瘤局部侵犯程度）
# 原始编码：1=submucosa, 2=muscle, 3=serosa, 4=contiguous
adtte$extent <- factor(adtte$extent,
                       levels = c(1, 2, 3, 4),
                       labels = c("Submucosa", "Muscle",
                                  "Serosa", "Contiguous"))

# 派生变量5：sex转因子
adtte$sex <- factor(adtte$sex,
                    levels = c(0, 1),
                    labels = c("Female", "Male"))


# ---- 3. 数据核查（Data Verification）-----------------------
# 对应SAP §9：Reproducibility
# 这一步必须跑，确认人群构造结果符合预期

cat("\n=== ITT人群核查 ===\n")
cat("ITT总N（预期619）:", nrow(adtte), "\n")

cat("\n治疗组分布:\n")
print(table(adtte$rx, dnn = "治疗组"))

cat("\n事件/删失分布:\n")
print(table(adtte$rx, adtte$status,
            dnn = c("治疗组", "status(1=事件,0=删失)")))

cat("\n分层因子node4分布:\n")
print(table(adtte$rx, adtte$node4,
            dnn = c("治疗组", "node4")))

cat("\n=== 关键协变量缺失情况 ===\n")
miss_vars <- c("nodes", "differ", "extent")
for (v in miss_vars) {
  n_miss <- sum(is.na(adtte[[v]]))
  cat(sprintf("  %-8s: NA = %d 人 (%.1f%%)\n",
              v, n_miss, 100 * n_miss / nrow(adtte)))
}

n_complete <- sum(complete.cases(adtte[, miss_vars]))
n_incomplete <- nrow(adtte) - n_complete
cat("\n任一协变量缺失（将被Cox自动剔除）:", n_incomplete, "人\n")
cat("完整案例（PP-like人群，预期594）:", n_complete, "人\n")


# ---- 4. 构造PP-like人群（用于SA-2）-------------------------
# 对应SAP §4：complete-case subset
# 显式标记，而非依赖Cox自动剔除

adtte_pp <- adtte[complete.cases(adtte[, miss_vars]), ]

cat("\n=== PP-like人群核查 ===\n")
cat("PP-like 总N（预期594）:", nrow(adtte_pp), "\n")
cat("PP-like 治疗组分布:\n")
print(table(adtte_pp$rx, dnn = "治疗组"))


# ---- 5. 保存分析数据集 --------------------------------------
# 在真实项目里，分析数据集要锁定（lock）后才能使用
# 这里我们保存为.RData，后续脚本直接load

save(adtte, adtte_pp,
     file = "adtte.RData")

cat("\n数据集已保存：adtte.RData\n")
cat("  adtte   : ITT人群，N =", nrow(adtte), "\n")
cat("  adtte_pp: PP-like人群，N =", nrow(adtte_pp), "\n")

cat("\n=== 3-A 完成 ===\n")

# ============================================================
# 项目：NCCTG结肠癌辅助化疗试验 生存分析复盘
# 脚本：3-B Table 1 基线特征表
# 对应SAP章节：§5 报告内容 / TFL: Table 1
# 输入：adtte.RData（由3-A生成）
# ============================================================


# ---- 1. 辅助函数：格式化输出 --------------------------------

# 连续变量：返回 "mean ± SD" 字符串
fmt_cont <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.1f ± %.1f", mean(x), sd(x))
}

# 分类变量某一水平：返回 "n (xx.x%)" 字符串
# 百分比分母用非缺失样本数（这是标准做法）
fmt_cat <- function(x, level) {
  x_nonmiss <- x[!is.na(x)]
  n <- sum(x_nonmiss == level)
  pct <- 100 * n / length(x_nonmiss)
  sprintf("%d (%.1f%%)", n, pct)
}

# 缺失数：返回 "n" 或在无缺失时返回 "0"
fmt_miss <- function(x) {
  as.character(sum(is.na(x)))
}


# ---- 2. 按治疗臂拆分数据 ------------------------------------
obs   <- adtte[adtte$rx == "Obs", ]
combo <- adtte[adtte$rx == "Lev+5FU", ]
all   <- adtte


# ---- 3. 逐行构造Table 1 -------------------------------------
# 每一行是一个 c(变量名, Obs列, Lev+5FU列, Overall列)

cat("============================================================\n")
cat("  Table 1. Baseline Characteristics by Treatment Arm\n")
cat("  (ITT Population, N = 619)\n")
cat("============================================================\n")

# 表头
hdr <- sprintf("%-32s %-18s %-18s %-18s",
               "Characteristic",
               sprintf("Obs (N=%d)", nrow(obs)),
               sprintf("Lev+5FU (N=%d)", nrow(combo)),
               sprintf("Overall (N=%d)", nrow(all)))
cat(hdr, "\n")
cat(strrep("-", 88), "\n")

# 行输出函数
prow <- function(label, c_obs, c_combo, c_all) {
  cat(sprintf("%-32s %-18s %-18s %-18s\n",
              label, c_obs, c_combo, c_all))
}

# --- Age（连续变量）---
prow("Age (years), mean ± SD",
     fmt_cont(obs$age), fmt_cont(combo$age), fmt_cont(all$age))

# --- Sex（分类）---
cat("Sex, n (%)\n")
prow("  Male",
     fmt_cat(obs$sex,"Male"), fmt_cat(combo$sex,"Male"),
     fmt_cat(all$sex,"Male"))
prow("  Female",
     fmt_cat(obs$sex,"Female"), fmt_cat(combo$sex,"Female"),
     fmt_cat(all$sex,"Female"))

# --- node4（分层因子）---
cat("Positive nodes, n (%)\n")
prow("  <=4 nodes",
     fmt_cat(obs$node4,"<=4 nodes"), fmt_cat(combo$node4,"<=4 nodes"),
     fmt_cat(all$node4,"<=4 nodes"))
prow("  >4 nodes",
     fmt_cat(obs$node4,">4 nodes"), fmt_cat(combo$node4,">4 nodes"),
     fmt_cat(all$node4,">4 nodes"))

# --- differ（肿瘤分化，有缺失）---
cat("Differentiation, n (%)\n")
for (lv in c("Well","Moderate","Poor")) {
  prow(paste0("  ", lv),
       fmt_cat(obs$differ,lv), fmt_cat(combo$differ,lv),
       fmt_cat(all$differ,lv))
}
prow("  Missing",
     fmt_miss(obs$differ), fmt_miss(combo$differ),
     fmt_miss(all$differ))

# --- extent（肿瘤侵犯程度）---
cat("Extent of invasion, n (%)\n")
for (lv in c("Submucosa","Muscle","Serosa","Contiguous")) {
  prow(paste0("  ", lv),
       fmt_cat(obs$extent,lv), fmt_cat(combo$extent,lv),
       fmt_cat(all$extent,lv))
}

# --- obstruct（结肠梗阻，0/1变量，作为附加基线特征）---
cat("Colon obstruction, n (%)\n")
prow("  Yes",
     fmt_cat(obs$obstruct,1), fmt_cat(combo$obstruct,1),
     fmt_cat(all$obstruct,1))

# --- 事件数（不是基线特征，但放表末供参考）---
cat(strrep("-", 88), "\n")
prow("Recurrence/death events, n (%)",
     fmt_cat(obs$status,1), fmt_cat(combo$status,1),
     fmt_cat(all$status,1))

cat("============================================================\n")
cat("注：百分比分母为非缺失样本数。\n")
cat("    按CONSORT规范，随机化试验Table 1不报组间检验p值。\n")
cat("============================================================\n")

cat("\n=== 3-B 完成 ===\n")
# ============================================================
# 项目：NCCTG结肠癌辅助化疗试验 生存分析复盘
# 脚本：3-C 主分析
#       KM曲线 + log-rank检验 + stratified Cox模型
# 对应SAP章节：§5.1 Primary Analysis
# 输入：adtte.RData（由3-A生成）
# 输出：Figure 1（KM曲线）、Table 2（主分析结果）
# ============================================================

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
# 输出到PDF（正式报告用）
pdf("Figure1_KM_curve.pdf", width = 8, height = 6)

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

dev.off()
cat("\nFigure 1 已保存：Figure1_KM_curve.pdf\n")

cat("\n=== 3-C 完成 ===\n")

# ============================================================
# 项目：NCCTG结肠癌辅助化疗试验 生存分析复盘
# 脚本：3-D 比例风险（PH）假设诊断
#       cox.zph检验 + Schoenfeld残差图
# 对应SAP章节：§5.1 Primary Analysis（模型假设核查）
# 输入：adtte.RData（由3-A生成）
# 输出：Figure 3（Schoenfeld残差图）
# ============================================================

# 重建3-C的主分析Cox模型
# 注意：node4是分层变量，分层变量本身不需要满足PH，
#       故cox.zph不检验node4，只检验 rx/age/sex/extent
cox_primary <- coxph(
  surv_obj ~ rx + strata(node4) + age + sex + extent,
  data = adtte,
  ties = "efron"
)

cat("主分析模型已重建，N =", cox_primary$n, "\n\n")


# ============================================================
# 第一部分：cox.zph() —— PH假设的统计检验
# ============================================================

# ---- 1. 运行cox.zph ----------------------------------------
# cox.zph对每个协变量检验"标准化Schoenfeld残差 ~ 时间"的斜率
#   H0: 斜率=0，即PH成立
#   p小 → 拒绝PH → 该变量违反比例风险假设
# transform参数：时间轴怎么变换后再做斜率检验
#   默认"km"——按KM估计变换，是常用稳健选择
ph_test <- cox.zph(cox_primary, transform = "km")

cat("=== cox.zph() PH假设检验结果 ===\n")
print(ph_test)

cat("\n--- 解读规则 ---\n")
cat("每行一个变量 + 最后GLOBAL总检验\n")
cat("p >= 0.05：无证据违反PH（假设大致成立）\n")
cat("p <  0.05：有证据违反PH（该变量效应随时间变化）\n")

# ---- 2. 逐变量判读 -----------------------------------------
cat("\n=== 逐变量判读 ===\n")
ph_tab <- ph_test$table
for (v in rownames(ph_tab)) {
  p_v <- ph_tab[v, "p"]
  verdict <- if (p_v >= 0.05) "PH大致成立" else "PH可能被违反"
  cat(sprintf("  %-10s: p = %.4f  -> %s\n", v, p_v, verdict))
}


# ============================================================
# 第二部分：Schoenfeld残差图 —— PH假设的图形诊断
# ============================================================

# ---- 3. 绘制Schoenfeld残差图 -------------------------------
# 每个变量一个子图：标准化Schoenfeld残差 vs 时间
#   - 散点：每个事件对应一个残差
#   - 实线：残差的平滑曲线，近似 beta(t) 的形状
#   - 虚线：平滑曲线的置信带
# 判读：
#   平滑线水平、置信带能包住一条水平线 → PH成立
#   平滑线明显倾斜/弯曲           → PH被违反
#
# ⚠️ "Schoenfeld残差图是PH诊断标准做法"为方法学共识，
#    建议核对 Therneau & Grambsch《Modeling Survival Data》

# 直接绘到当前图形设备：本地运行时即 R/RStudio 图形窗口，
# 图会即时弹出，便于自己观察、形成判断（不写pdf()堵住窗口）

# 子图布局：rx/age/sex/extent 共4个变量 → 2x2
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

# plot.cox.zph 会为每个变量画一张Schoenfeld残差图
plot(ph_test,
     resid = TRUE,    # 显示残差散点
     se    = TRUE,    # 显示置信带
     col   = c("firebrick", "steelblue"))  # 平滑线/置信带颜色

# 恢复单图布局
par(mfrow = c(1, 1))

cat("\nSchoenfeld残差图已绘制到图形窗口。\n")
cat("请自行观察4条平滑线是否水平、置信带能否包住水平线。\n")

# ---- 可选：如需存文件，解开下面注释 -------------------------
# pdf("Figure3_Schoenfeld.pdf", width = 9, height = 7)
# par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
# plot(ph_test, resid = TRUE, se = TRUE,
#      col = c("firebrick", "steelblue"))
# dev.off()


# ============================================================
# 第三部分：结果汇总与结论
# ============================================================

# ---- 4. 生成结论文字 ---------------------------------------
global_p <- ph_tab["GLOBAL", "p"]

cat("\n=== PH诊断结论 ===\n")
cat(sprintf("GLOBAL检验 p = %.4f\n", global_p))

if (global_p >= 0.05) {
  cat("结论：GLOBAL检验未拒绝PH假设。\n")
  cat("      主分析Cox模型的比例风险假设大致成立，\n")
  cat("      HR=0.59可解释为整个随访期内恒定的风险比。\n")
} else {
  cat("结论：GLOBAL检验拒绝PH假设，需逐变量排查。\n")
  cat("      若违反来自治疗变量rx，主分析HR应解释为\n")
  cat("      时间平均效应，并启用预先指定的RMST敏感性分析。\n")
}

# ---- 5. 特别检查治疗变量rx ---------------------------------
# rx是研究变量，它的PH是否成立最关键
p_rx <- ph_tab["rx", "p"]
cat("\n=== 治疗变量rx的PH检验（最关键）===\n")
cat(sprintf("rx 的 PH检验 p = %.4f\n", p_rx))
if (p_rx >= 0.05) {
  cat("治疗效应的比例风险假设成立——\n")
  cat("HR=0.59是对治疗效应的恰当总结。\n")
} else {
  cat("治疗效应可能随时间变化——\n")
  cat("需在CSR中说明，并参考RMST敏感性分析。\n")
}

cat("\n=== 3-D 完成 ===\n")

