# ============================================================
# 项目：NCCTG结肠癌辅助化疗试验 生存分析复盘
# 脚本：3-D 比例风险（PH）假设诊断
#       cox.zph检验 + Schoenfeld残差图
# 对应SAP章节：§5.1 Primary Analysis（模型假设核查）
# 输入：adtte.RData（由3-A生成）
# 输出：Figure 3（Schoenfeld残差图）
# ============================================================

rm(list = ls())
library(survival)

# ---- 0. 载入数据并重建主分析模型 ----------------------------
load("adtte.RData")
surv_obj <- Surv(adtte$time_months, adtte$status)

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
  cat("      HR=0.62可解释为整个随访期内恒定的风险比。\n")
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
  cat("HR=0.62是对治疗效应的恰当总结。\n")
} else {
  cat("治疗效应可能随时间变化——\n")
  cat("需在CSR中说明，并参考RMST敏感性分析。\n")
}

cat("\n=== 3-D 完成 ===\n")
