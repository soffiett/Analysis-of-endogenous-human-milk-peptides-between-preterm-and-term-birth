Pep = read.csv ("Peptide_Profile.csv",header = FALSE,stringsAsFactors=FALSE)
Pep = Pep[-c(1, 222,263),]
vars = data.frame(cbind (maturation, time,mum_id))
Pepmat =  Pep [,6:65]

t1 <- apply(Pepmat[,vars$maturation=="1"]==1,1,sum)
t2 <- apply(Pepmat[,vars$maturation=="2"]==1,1,sum)

peptides.pretermonly = Pep[t2==32,5]
peptides.termonly = Pep[t1==28,5]


ind<- (t1 <= 21)&(t2 <= 24)

Pepmat.clean  = Pepmat[t1<=21&t2<=24,] #omit the peptide that are absent in 75% of sample
Pep.clean = Pep [t1<=21&t2<=24,]
Pepmat.log = log (Pepmat.clean)

library (lme4)
#Main effect
lms.main <- function(mat,vars)
{
  m <- nrow(mat)
  pvmat <- matrix(rep(0,m),ncol=1)
  for (i in 1:m)
  {
    new.dat <- data.frame(vars,t(mat[i,]))
    names(new.dat)[4] <- "y"
    aobj <- anova(lmer(y~maturation:time+(1|mum_id), data = new.dat, REML = FALSE), 
                  lmer(y~maturation+maturation:time+(1|mum_id), data = new.dat,REML = FALSE))
    pvs <- aobj$Pr[2]
    pvmat[i] <- pvs
  }
  return(pvmat)
}


p.main = lms.main (Pepmat.log, vars)
p.main = p.adjust (p.main, "fdr")

require("MASS")
truehist(p.main)
table(p.main <0.05)
Pep.result = data.frame (identity = Pep.clean [p.main<0.05,5], p_main=p.main [p.main<0.05])
colnames(Pep.result) = c("Peptide_identity","P-value_of_Likelihood_Test_For_Main")

#Interaction effect
lms.int <- function(mat,vars)
{
  m <- nrow(mat)
  pvmat <- matrix(rep(0,m),ncol=1)
  for (i in 1:m)
  {
    new.dat <- data.frame(vars,t(mat[i,]))
    names(new.dat)[4] <- "y"
    aobj <- anova(lmer(y~maturation+ (1|mum_id), data = new.dat,REML = FALSE),
                  lmer(y~maturation+ maturation: time+(1|mum_id), data = new.dat,REML = FALSE))
    pvs <- aobj$Pr[2]
    pvmat[i] <- pvs
  }
  return(pvmat)
}

p.int = lms.int (Pepmat.log, vars)
p.int = p.adjust (p.int, "fdr")

table(p.int <0.05)

Pep.result2 = data.frame (Pep.clean [p.int<0.05,5], p.int [p.int<0.05])
colnames(Pep.result2) = c("Peptide_identity","P-value_of_Likelihood_Test_For_Interaction")

#Separate term vs. preterm
pre.mat = Pepmat.log[,1:28]
term.mat = Pepmat.log[,29:60]

lms.term <- function(mat,vars)
{
  m <- nrow(mat)
  pvmat <- matrix(rep(0,m),ncol=1)
  for (i in 1:m)
  {
    new.dat <- data.frame(vars,t(mat[i,]))
    names(new.dat)[4] <- "y"
    aobj <- anova(lmer(y~time+(1|mum_id), data = new.dat, REML=FALSE), 
                         lmer(y~1+ (1|mum_id), data = new.dat,REML=FALSE))
    pvs <- aobj$Pr[2]
    pvmat[i] <- pvs
  }
  return(pvmat)
}
p.pre = lms.term (pre.mat, vars[1:28,])
p.pre = p.adjust (p.pre, "fdr") #

p.term = lms.term (term.mat, vars[29:60,])
p.term = p.adjust (p.term, "fdr") #

preterm.result = data.frame (Pep.clean [p.pre<0.05,5], p.pre[p.pre<0.05])
term.result = data.frame (Pep.clean [p.term<0.05,5], p.term[p.term<0.05])
colnames(preterm.result) = c("Peptide_identity","P-value_of_F_test")
colnames(term.result) = c("Peptide_identity","P-value_of_F_test")

library("xlsx")
wb = createWorkbook()          # create blank workbook
sheet1 = createSheet(wb, sheetName="Significant Peptide in Main Effect") # create different sheets
addDataFrame(Pep.result, sheet1)  # add data to the sheets

sheet2 = createSheet(wb, sheetName="Significant Peptide in Interaction Effect") # create different sheets
addDataFrame(Pep.result2, sheet2)  # add data to the sheets

sheet3 = createSheet(wb, sheetName="Significant Peptides in Preterm Samples") # create different sheets
addDataFrame(preterm.result, sheet3)  # add data to the sheets

sheet4 = createSheet(wb, sheetName="Significant Peptides in Term Samples") # create different sheets
addDataFrame(term.result, sheet4)  # add data to the sheets


sheet5 = createSheet(wb, sheetName="Peptides only shown in preterm samples") # create different sheets
addDataFrame(peptides.pretermonly, sheet5)  # add data to the sheets

sheet6 = createSheet(wb, sheetName="Peptides only shown in term samples") # create different sheets
addDataFrame(peptides.termonly, sheet6)  # add data to the sheets

saveWorkbook(wb, "Statistical Analysis of Peptides-102114.xlsx")  # write the file with multiple sheets


#########################################
lms.main2 <- function(mat,vars)
{
  m <- nrow(mat)
  pvmat <- matrix(rep(0,2*m),ncol=2)
  for (i in 1:m)
  {
    new.dat <- data.frame(vars,t(mat[i,]))
    names(new.dat)[3] <- "y"
    aobj <- anova(lm(y~maturation+time, data = new.dat))
    pvs <- aobj$Pr[-3]
    pvmat[i,] <- pvs
  }
  return(pvmat)
}

p.main2 = lms.main2 (Pepmat.log, vars)
p.main2 = apply (p.main2, 2, function (x) (p.adjust(x, "fdr")))

apply (p.main2, 2, function (x) sum(x<=0.05))

#check omit peptide

Pepmat.omit = Pepmat[onecounts>=45,]
#heatmap
Pep.heat = Pep.clean [p.main<=0.05,]
heat.identity = Pep.heat[,5]

heat.mean= matrix (NA, ncol = 2, nrow = nrow (Pep.heat))
colnames(heat.mean) = c("Preterm","Term")

for (i in 1: nrow (Pep.heat)) {
  
  heat.mean[i,]= tapply(t(Pep.heat [i,-c(1:5)]), maturation, mean)
}



Pep.heat.mat = data.matrix(heat.mean)
rownames(Pep.heat.mat) = heat.identity

heatmap(Pep.heat.mat, Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))



