peptide = read.csv ("peptide.csv", header = FALSE)
p.amount = peptide[,4]
maturation = peptide[,1]
time = peptide[,2]
mum_id =as.factor( peptide[,3])

with(peptide, interaction.plot(time, maturation, p.amount,
    ylab = "mean of p.amount", xlab = "time", trace.label = "maturation"))

gobj = glm(p.amount~maturation+time,family=quasipoisson)
drop1(gobj,test="Chisq")

anova(gobj, test = "Chisq")

gobj2 = glm(p.amount~maturation*time+Error(mum_id),family=quasipoisson)
anova(gobj2, test = "Chisq")
summary(aov(p.amount~maturation*time+Error(mum_id))
        
 
#=========================using the mix effect model============================#
fm1 = glmer (p.amount ~ maturation + time + maturation: time + 
                 (1|mum_id), family = poisson)  #4.28e-10 for maturation,  7.57e-05 for time D, 7.56e-13 *** for interaction
#notice that the p value from output of glmer can only be used as general guide
#to test your model, you need to do likelihood ratio test for the real p value
fm2 = glmer (p.amount ~ maturation + time  + 
                  (1|mum_id), family = poisson)   
fm3 = glmer (p.amount ~  time + maturation: time + 
                  (1|mum_id), family = poisson)    
fm4 = glmer (p.amount ~ maturation +  maturation: time + 
                (1|mum_id), family = poisson)     
#conduct likelihood test
anova(fm2, fm) #<2.2*e-16
anova(fm3, fm) #<2.2*e-16
anova(fm4, fm) #<2.2*e-16
#conclusion: maturation, last lactation stage, interaction between
#  maturation and last lactation stage

#For peptide intensity
Pep = read.csv ("Peptide_Profile.csv",header = FALSE,stringsAsFactors=FALSE)
p.int = apply (Pep[,-c(1:5)], 2, mean)

fm5 = lmer (p.int ~ maturation + time + maturation: time + (1|mum_id)) 
fm6 = lmer (p.int ~ maturation + time +  (1|mum_id)) #test interaction fact
anova(fm5,fm6) #2.863e-05 *** for interaction
fm7 = lmer (p.int ~  time + maturation: time + (1|mum_id)) #test maturation
anova(fm5,fm7) #<2.2*e-16
fm8 = lmer (p.int ~ maturation  + maturation: time + (1|mum_id)) #test effect of time
anova (fm5,fm8) 

#For pep concentration
p.conc = read.csv ("BCA_con.csv", header = FALSE, stringsAsFactors = FALSE)
p.conc = as.numeric (unlist(p.conc))
fm9 = lmer (p.conc ~ maturation[-40] + time[-40] + maturation[-40] : time[-40] + (1|mum_id[-40])) 
fm10 = lmer (p.conc ~  time[-40] + maturation[-40]: time[-40] + (1|mum_id[-40]))
fm11 = lmer (p.conc ~  maturation[-40] + maturation[-40] : time[-40]+(1|mum_id[-40]))
fm12 = lmer (p.conc ~ maturation[-40] + time[-40] + (1|mum_id[-40])) 

#For Plasmin Activity
p.plas = read.csv ("plasmin.csv", header = FALSE, stringsAsFactors = FALSE)
p.plas = as.numeric (unlist(p.plas))
fm13 = lmer (p.plas ~  maturation + time  + (1|mum_id), REML=FALSE) 
fm14 = lmer (p.plas ~  time + maturation: time + (1|mum_id),REML=FALSE)
fm15 = lmer (p.plas ~   maturation + (1|mum_id),REML=FALSE)
fm16 = lmer (p.plas ~ maturation + time + (1|mum_id),REML=FALSE) 


#for protein analysis
by_protein = read.csv ("by_protein.csv", header = TRUE)

fm17 = lmer (by_protein[,1] ~ maturation + time+maturation:time + (1|mum_id))#, REML=FALSE) 
fm18 = lmer (by_protein[,1] ~  maturation + time:maturation +  (1|mum_id))#,REML=FALSE)

fm19 = lmer (by_protein[,2] ~  maturation + time +(1|mum_id))
fm20 = lmer (by_protein[,2] ~ maturation +  (1|mum_id)) 

fm21 = lmer (by_protein[,3] ~   time + maturation:time+ (1|mum_id))
fm22 = lmer (by_protein[,3] ~ 1 +  maturation:time  +(1|mum_id)) 

fm23 = lmer (by_protein[,4] ~  maturation + time  + (1|mum_id) , REML= FALSE)
fm24 = lmer (by_protein[,4] ~ maturation  +  (1|mum_id), REML = FALSE) 


fm25 = lmer (by_protein[,5] ~  maturation + time  +(1|mum_id))
fm26 = lmer (by_protein[,5] ~ maturation +  (1|mum_id)) 


p.values.lmer <- function(x) {
  summary.model <- summary(x)
  data.lmer <- data.frame(model.matrix(x))
  names(data.lmer) <- names(fixef(x))
  names(data.lmer) <- gsub(pattern=":", x=names(data.lmer), replacement=".", fixed=T)
  names(data.lmer) <- ifelse(names(data.lmer)=="(Intercept)", "Intercept", names(data.lmer))
  string.call <- strsplit(x=as.character(x@call), split=" + (", fixed=T)
  var.dep <- unlist(strsplit(x=unlist(string.call)[2], " ~ ", fixed=T))[1]
  vars.fixef <- names(data.lmer)
  formula.ranef <- paste("+ (", string.call[[2]][-1], sep="")
  formula.ranef <- paste(formula.ranef, collapse=" ")
  formula.full <- as.formula(paste(var.dep, "~ -1 +", paste(vars.fixef, collapse=" + "), 
                                   formula.ranef))
  data.ranef <- data.frame(x@frame[, 
                                   which(names(x@frame) %in% names(ranef(x)))])
  names(data.ranef) <- names(ranef(x))
  data.lmer <- data.frame(x@frame[, 1], data.lmer, data.ranef)
  names(data.lmer)[1] <- var.dep
  out.full <- lmer(formula.full, data=data.lmer, REML=F)
  p.value.LRT <- vector(length=length(vars.fixef))
  for(i in 1:length(vars.fixef)) {
    formula.reduced <- as.formula(paste(var.dep, "~ -1 +", paste(vars.fixef[-i], 
                                                                 collapse=" + "), formula.ranef))
    out.reduced <- lmer(formula.reduced, data=data.lmer, REML=F)
    print(paste("Reduced by:", vars.fixef[i]))
    print(out.LRT <- data.frame(anova(out.full, out.reduced)))
    p.value.LRT[i] <- round(out.LRT[2, 7], 3)
  }
  summary.model@coefs <- cbind(summary.model@coefs, p.value.LRT)
  summary.model@methTitle <- c("\n", summary.model@methTitle, 
                               "\n(p-values from comparing nested models fit by maximum likelihood)")
  print(summary.model)
}

library(lme4)
library(SASmixed)
lmer.out <- lmer(strength ~ Program * Time + (Time|Subj), data=Weights)
p.values.lmer(lmer.out)

coefs <- data.frame(coef(summary(trial2)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs

require(pbkrtest)
df.KR <- get_ddf_Lb(trial2, fixef(trial2))
# get p-values from the t-distribution using the t-values and approximated
# degrees of freedom
coefs$p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
coefs

fmp1 = lmer (by_protein[,1] ~ maturation + time +  maturation:time +  (1|mum_id))
df.KR1 <- get_ddf_Lb(fmp1, fixef(fmp1))
coefs1 <- data.frame(coef(summary(fmp1)))
# get p-values from the t-distribution using the t-values and approximated
# degrees of freedom
coefs1$p.KR <- 2 * (1 - pt(abs(coefs1$t.value), df.KR1))
p.adjust (coefs1$p.KR,"fdr")

fmp2 = lmer (by_protein[,2] ~ maturation + time +  maturation:time +  (1|mum_id))
df.KR2 <- get_ddf_Lb(fmp2, fixef(fmp2))
coefs2 <- data.frame(coef(summary(fmp2)))
# get p-values from the t-distribution using the t-values and approximated
# degrees of freedom
coefs2$p.KR <- 2 * (1 - pt(abs(coefs2$t.value), df.KR2))
coefs2

fmp3 = lmer (by_protein[,3] ~ maturation + time +  maturation:time +  (1|mum_id))
df.KR3 <- get_ddf_Lb(fmp3, fixef(fmp3))
coefs3 <- data.frame(coef(summary(fmp3)))
# get p-values from the t-distribution using the t-values and approximated
# degrees of freedom
coefs3$p.KR <- 2 * (1 - pt(abs(coefs3$t.value), df.KR3))
coefs3


fmp4 = lmer (by_protein[,4] ~ maturation + time +  maturation:time +  (1|mum_id))
df.KR4 <- get_ddf_Lb(fmp4, fixef(fmp4))
coefs4 <- data.frame(coef(summary(fmp4)))
# get p-values from the t-distribution using the t-values and approximated
# degrees of freedom
coefs4$p.KR <- 2 * (1 - pt(abs(coefs4$t.value), df.KR4))
coefs4

fmp5 = lmer (by_protein[,5] ~ maturation + time +  maturation:time +  (1|mum_id))
df.KR5 <- get_ddf_Lb(fmp5, fixef(fmp5))
coefs5 <- data.frame(coef(summary(fmp5)))
# get p-values from the t-distribution using the t-values and approximated
# degrees of freedom
coefs5$p.KR <- 2 * (1 - pt(abs(coefs5$t.value), df.KR5))
coefs5

fmp6 = lmer (p.plas ~ maturation + time +  maturation:time +  (1|mum_id))
df.KR6 <- get_ddf_Lb(fmp6, fixef(fmp6))
coefs6 <- data.frame(coef(summary(fmp6)))
# get p-values from the t-distribution using the t-values and approximated
# degrees of freedom
coefs6$p.KR <- 2 * (1 - pt(abs(coefs6$t.value), df.KR6))
coefs6