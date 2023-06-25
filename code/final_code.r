ir_Dsec_TW_F = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\ir_Dsec_TW_F.csv', header=TRUE, sep=",")
ir_Dsec_TW_M = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\ir_Dsec_TW_M.csv', header=TRUE, sep=",")
obp_Dsec_TW_F = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\obp_Dsec_TW_F.csv', header=TRUE, sep=",")
obp_Dsec_TW_M = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\obp_Dsec_TW_M.csv', header=TRUE, sep=",")
or_Dsec_TW_F = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\or_Dsec_TW_F.csv', header=TRUE, sep=",")
or_Dsec_TW_M = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\or_Dsec_TW_M.csv', header=TRUE, sep=",")
ir_Dsec_JP_F = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\ir_Dsec_JP_F.csv', header=TRUE, sep=",")
ir_Dsec_JP_M = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\ir_Dsec_JP_M.csv', header=TRUE, sep=",")
obp_Dsec_JP_F = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\obp_Dsec_JP_F.csv', header=TRUE, sep=",")
obp_Dsec_JP_M = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\obp_Dsec_JP_M.csv', header=TRUE, sep=",")
or_Dsec_JP_F = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\or_Dsec_JP_F.csv', header=TRUE, sep=",")
or_Dsec_JP_M = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\or_Dsec_JP_M.csv', header=TRUE, sep=",")
ir_Dsim_JP_F = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\ir_Dsim_JP_F.csv', header=TRUE, sep=",")
ir_Dsim_JP_M = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\ir_Dsim_JP_M.csv', header=TRUE, sep=",")
obp_Dsim_JP_F = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\obp_Dsim_JP_F.csv', header=TRUE, sep=",")
obp_Dsim_JP_M = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\obp_Dsim_JP_M.csv', header=TRUE, sep=",")
or_Dsim_JP_F = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\or_Dsim_JP_F.csv', header=TRUE, sep=",")
or_Dsim_JP_M = read.table('C:\\Users\\ASUS\\Desktop\\bioinfo\\final_project\\no Orco no Dmel\\or_Dsim_JP_M.csv', header=TRUE, sep=",")


heat_ir <- data.frame(Dsec_F_TW=c(ir_Dsec_TW_F[,5]),
                      Dsec_F_JP=c(ir_Dsec_JP_F[,5]),
                      Dsec_M_TW=c(ir_Dsec_TW_M[,5]),
                      Dsec_M_JP=c(ir_Dsec_JP_M[,5]),
                      Dsim_F_JP=c(ir_Dsim_JP_F[,5]),
                      Dsim_M_JP=c(ir_Dsim_JP_M[,5]))
rownames(heat_ir) = ir_Dsec_TW_F[,1]
heat_obp <- data.frame(Dsec_F_TW=c(obp_Dsec_TW_F[,5]),
                      Dsec_F_JP=c(obp_Dsec_JP_F[,5]),
                      Dsec_M_TW=c(obp_Dsec_TW_M[,5]),
                      Dsec_M_JP=c(obp_Dsec_JP_M[,5]),
                      Dsim_F_JP=c(obp_Dsim_JP_F[,5]),
                      Dsim_M_JP=c(obp_Dsim_JP_M[,5]))
rownames(heat_obp) = obp_Dsec_TW_F[,1]
heat_or <- data.frame(Dsec_F_TW=c(or_Dsec_TW_F[,5]),
                      Dsec_F_JP=c(or_Dsec_JP_F[,5]),
                      Dsec_M_TW=c(or_Dsec_TW_M[,5]),
                      Dsec_M_JP=c(or_Dsec_JP_M[,5]),
                      Dsim_F_JP=c(or_Dsim_JP_F[,5]),
                      Dsim_M_JP=c(or_Dsim_JP_M[,5]))
rownames(heat_or) = or_Dsec_TW_F[,1]

heat_ir_scale <- scale(heat_ir)
heat_obp_scale <- scale(heat_obp)
heat_or_scale <- scale(heat_or)
#heatmap(x, scale='none')


library('pvclust')
ir_bp <- pvclust(heat_ir_scale, nboot=100)
obp_bp <- pvclust(heat_obp_scale, nboot=100)
or_bp <- pvclust(heat_or_scale, nboot=100)
#plot(x)

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
expMatrix_ir <- heat_ir
expMatrix_obp <- heat_obp
expMatrix_or <- heat_or

tpms_ir <- apply(expMatrix_ir, 2, fpkmToTpm)
tpms_obp <- apply(expMatrix_obp, 2, fpkmToTpm)
tpms_or <- apply(expMatrix_or, 2, fpkmToTpm)


group_list = c(rep('Normal', 3), rep('Tumor', 3))
## 強制限定順序
group_list <- factor(group_list, levels = c("Normal", "Tumor"), ordered = F)
#表達矩陣數據校正
exprSet_ir <- tpms_ir
exprSet_obp <- tpms_obp
exprSet_or <- tpms_or

library(limma) 
exprSet_ir = normalizeBetweenArrays(exprSet_ir)
exprSet_obp = normalizeBetweenArrays(exprSet_obp)
exprSet_or = normalizeBetweenArrays(exprSet_or)

#判斷數據是否需要轉換
exprSet_ir <- log2(exprSet_ir+1)
exprSet_obp <- log2(exprSet_obp+1)
exprSet_or <- log2(exprSet_or+1)
#差異分析：
dat_ir <- exprSet_ir
dat_obp <- exprSet_obp
dat_or <- exprSet_or
design = model.matrix(~factor(group_list))
fit_ir = lmFit(dat_ir, design)
fit_obp = lmFit(dat_obp, design)
fit_or = lmFit(dat_or, design)
fit_ir = eBayes(fit_ir)
fit_obp = eBayes(fit_obp)
fit_or = eBayes(fit_or)
options(digits = 4)

deg_ir = topTable(fit_ir, coef=2, adjust='BH', number = Inf)
deg_obp = topTable(fit_obp, coef=2, adjust='BH', number = Inf)
deg_or = topTable(fit_or, coef=2, adjust='BH', number = Inf)