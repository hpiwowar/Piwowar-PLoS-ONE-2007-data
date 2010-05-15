
# Read in data
setwd("e:\\research\\PLoS ONE submission\\Revised\\")
dat = read.csv("PiwowarData.csv", header=TRUE, row.names=1)

dim(dat)
rownames(dat)
colnames(dat)

# Calculate a few fields which will be useful later
dat$cohortmonthsfromend =
  max(dat$Number.of.months.between.1.99.and.trial.publication) -
  dat$Number.of.months.between.1.99.and.trial.publication
dat$Number.of.cases.in.trial.gt25 = dat$Number.of.cases.in.trial > 25


# A quick summary of the data
print("Number of papers")
print("Data not available, Data available")
tapply(dat$Number.of.Citations.during.2004.2005 > 1,
       dat$Is.the.microarray.data.publicly.available,
       sum)

print("Number of citations")
print("Data not available, Data available")
tapply(dat$Number.of.Citations.during.2004.2005,
       dat$Is.the.microarray.data.publicly.available,
       sum)



## Table 1

# Note that the Fisher estimates of the odd's ratio are not exactly the same as a*d/b*c.  The paper actually reports the latter, but uses the fisher values for the confidence interval

impact.2x2 	= table(dat$Impact.factor.of.journal < 25, !dat$Is.the.microarray.data.publicly.available)
print(impact.2x2)
impact.fisher = fisher.test(impact.2x2)
print(round(impact.fisher$estimate, 1))
print(round(impact.fisher$conf.int, 1))

year.2x2 	= table(dat$Number.of.months.between.1.99.and.trial.publication > 24, !dat$Is.the.microarray.data.publicly.available)
print(year.2x2)
year.fisher = fisher.test(year.2x2)
print(round(year.fisher$estimate, 1))
print(round(year.fisher$conf.int, 1))

usauth.2x2 	= table(dat$Are.there.any.authors.from.the.US == 0, !dat$Is.the.microarray.data.publicly.available)
print(usauth.2x2)
usauth.fisher = fisher.test(usauth.2x2)
print(round(usauth.fisher$estimate, 1))
print(round(usauth.fisher$conf.int, 1))





### Table 2
## Primary analysis


# Some helper functions
calcCI.exp= function(res, param) {
  coefs = summary(res)$coeff
  coeff = coefs[param,]
  x = coeff[1]
  stderr = coeff[2]
  p = coeff[4]
  return(list(param = param,
              est = round(exp(x), 2), 
	      CI = c(round(exp(x - 1.96*stderr), 2),
                     round(exp(x + 1.96*stderr), 2)), 
  	      p = round(p, 3)))
}

calcCI.noexp= function(res, param) {
  coefs = summary(res)$coeff
  coeff = coefs[param,]
  x = coeff[1]
  stderr = coeff[2]
  p = coeff[4]
  return(list(param = param,
              est = round(x, 2), 
	      CI = c(round(x - 1.96*stderr, 2),
                     round(x + 1.96*stderr, 2)), 
  	      p = round(p, 3)))
}

all.results = function(res) {
  # give the results of the impact factor without exp because it is the
  # log impact factor, so interpretation is easier if kept in the log domain
  print(calcCI.noexp(res, "lnimpact"))
  print(calcCI.exp(res, "Are.there.any.authors.from.the.US"))
  print(calcCI.exp(res, "Number.of.months.between.1.99.and.trial.publication"))
  print(calcCI.exp(res, "Is.the.microarray.data.publicly.available"))
}


# Take the log of the endpoints and impact factor
dat$lnimpact = log(dat$Impact.factor.of.journal)
dat$lncites0405 = log(dat$Number.of.Citations.during.2004.2005)
dat$lncites24months = log(dat$Number.of.Citations.in.first.24.months.after.publication)

# Define the lower-profile subset
which.subset = which((dat$Impact.factor.of.journal < 25) & (dat$Number.of.months.between.1.99.and.trial.publication > 24))
dat.subset = dat[which.subset,]

# A quick summary of the lower-profile subset
print("Number of papers in the lower-profile subset")
print("Data not available, Data available")
tapply(dat.subset$Number.of.Citations.during.2004.2005 > 1,
       dat.subset$Is.the.microarray.data.publicly.available,
       sum)

print("Number of citations in the lower profile-subset")
print("Data not available, Data available")
tapply(dat.subset$Number.of.Citations.during.2004.2005,
       dat.subset$Is.the.microarray.data.publicly.available,
       sum)


# Do the regressions
result.primary = lm(lncites0405
  ~ lnimpact + Number.of.months.between.1.99.and.trial.publication + Are.there.any.authors.from.the.US + Is.the.microarray.data.publicly.available, dat)
all.results(result.primary)

result.primary.24mo = lm(lncites24months
  ~ lnimpact + Number.of.months.between.1.99.and.trial.publication + Are.there.any.authors.from.the.US + Is.the.microarray.data.publicly.available, dat)
all.results(result.primary.24mo)

result.primary.subset = lm(lncites0405
  ~ lnimpact + Number.of.months.between.1.99.and.trial.publication + Are.there.any.authors.from.the.US + Is.the.microarray.data.publicly.available, dat.subset)
all.results(result.primary.subset)


# Table 3
# Exploratory results
# Articles, No. is confirmed below
# Citations, No. is therefore not confirmed; assumed to be correct based on prior computation
# coef is raised to 10 and called "fold increase" since covariates are binary
# note the new p-values

# Use subset that makes data available
which.Is.the.microarray.data.publicly.available = which(dat$Is.the.microarray.data.publicly.available == 1)
dat.da = dat[which.Is.the.microarray.data.publicly.available,]
n.Is.the.microarray.data.publicly.available = length(which.Is.the.microarray.data.publicly.available)

# Do the calculations
result.expl.n 			= lm(lncites0405
  ~ Number.of.cases.in.trial.gt25 + lnimpact + Are.there.any.authors.from.the.US + Number.of.months.between.1.99.and.trial.publication, 
  subset=which.Is.the.microarray.data.publicly.available, dat)
print(result.expl.n$call)
print(round(exp(summary(result.expl.n)$coeff[2,]), 2))
tab = table(dat.da$Number.of.cases.in.trial.gt25); tab; round(tab/sum(tab), 2)
tab = tapply(dat.da$Number.of.Citations.during.2004.2005, dat.da$Number.of.cases.in.trial.gt25, sum); tab; round(tab/sum(tab), 2)

result.expl.clinical		= lm(lncites0405
  ~ Trial.has.a.clinical.endpoint + lnimpact + Are.there.any.authors.from.the.US + Number.of.months.between.1.99.and.trial.publication, 
  subset=which.Is.the.microarray.data.publicly.available, dat)
print(result.expl.clinical$call)
print(round(exp(summary(result.expl.clinical)$coeff[2,]), 2))
tab = table(dat.da$Trial.has.a.clinical.endpoint); tab; round(tab/sum(tab), 2)
tab = tapply(dat.da$Number.of.Citations.during.2004.2005, dat.da$Trial.has.a.clinical.endpoint, sum); tab; round(tab/sum(tab), 2)

result.expl.affy		= lm(lncites0405
  ~ Uses.the.Affymetrix.microarray.platform + lnimpact + Are.there.any.authors.from.the.US + Number.of.months.between.1.99.and.trial.publication, 
  subset=which.Is.the.microarray.data.publicly.available, dat)
print(result.expl.affy$call)
print(round(exp(summary(result.expl.affy)$coeff[2,]), 2))
tab = table(dat.da$Uses.the.Affymetrix.microarray.platform); tab; round(tab/sum(tab), 2)
tab = tapply(dat.da$Number.of.Citations.during.2004.2005, dat.da$Uses.the.Affymetrix.microarray.platform, sum); tab; round(tab/sum(tab), 2)

result.expl.geo			= lm(lncites0405
  ~ In.the.GEO.database + lnimpact + Are.there.any.authors.from.the.US + Number.of.months.between.1.99.and.trial.publication, 
  subset=which.Is.the.microarray.data.publicly.available, dat)
print(result.expl.geo$call)
print(round(exp(summary(result.expl.geo)$coeff[2,]), 2))
tab = table(dat.da$In.the.GEO.database); tab; round(tab/sum(tab), 2)
tab = tapply(dat.da$Number.of.Citations.during.2004.2005, dat.da$In.the.GEO.database, sum); tab; round(tab/sum(tab), 2)

result.expl.smd			= lm(lncites0405
  ~ In.the.SMD.database + lnimpact + Are.there.any.authors.from.the.US + Number.of.months.between.1.99.and.trial.publication, 
  subset=which.Is.the.microarray.data.publicly.available, dat)
print(result.expl.smd$call)
print(round(exp(summary(result.expl.smd)$coeff[2,]), 2))
tab = table(dat.da$In.the.SMD.database); tab; round(tab/sum(tab), 2)
tab = tapply(dat.da$Number.of.Citations.during.2004.2005, dat.da$In.the.SMD.database, sum); tab; round(tab/sum(tab), 2)

result.expl.raw			= lm(lncites0405
  ~ Raw.data.such.as.CEL.files.are.available + lnimpact + Are.there.any.authors.from.the.US + Number.of.months.between.1.99.and.trial.publication,
  subset=which.Is.the.microarray.data.publicly.available, dat)
print(result.expl.raw$call)
print(round(exp(summary(result.expl.raw)$coeff[2,]), 2))
tab = table(dat.da$Raw.data.such.as.CEL.files.are.available); tab; round(tab/sum(tab), 2)
tab = tapply(dat.da$Number.of.Citations.during.2004.2005, dat.da$Raw.data.such.as.CEL.files.are.available, sum); tab; round(tab/sum(tab), 2)

result.expl.suppl		= lm(lncites0405
  ~ Publication.mentions.Supplemental.data + lnimpact + Are.there.any.authors.from.the.US + Number.of.months.between.1.99.and.trial.publication,
  subset=which.Is.the.microarray.data.publicly.available, dat)
print(result.expl.suppl$call)
print(round(exp(summary(result.expl.suppl)$coeff[2,]), 2))
tab = table(dat.da$Publication.mentions.Supplemental.data); tab; round(tab/sum(tab), 2)
tab = tapply(dat.da$Number.of.Citations.during.2004.2005, dat.da$Publication.mentions.Supplemental.data, sum); tab; round(tab/sum(tab), 2)

result.expl.oncomine		= lm(lncites0405
  ~ Publication.has.an.Oncomine.profile + lnimpact + Are.there.any.authors.from.the.US + Number.of.months.between.1.99.and.trial.publication,
  subset=which.Is.the.microarray.data.publicly.available, dat)
print(result.expl.oncomine$call)
print(round(exp(summary(result.expl.oncomine)$coeff[2,]), 2))
tab = table(dat.da$Publication.has.an.Oncomine.profile); tab; round(tab/sum(tab), 2)
tab = tapply(dat.da$Number.of.Citations.during.2004.2005, dat.da$Publication.has.an.Oncomine.profile, sum); tab; round(tab/sum(tab), 2)


### Figure 1

table(dat$Is.the.microarray.data.publicly.available)
boxplot(Number.of.Citations.during.2004.2005 ~ Is.the.microarray.data.publicly.available,
        data = dat,
        boxwex = 0.5, 
        names=c("Data Not Shared (n=44)", "Data Shared (n=41)"), 
        ylab = "Number of Citations in 2004-2005", outline=T, notch=F, log="y")
dev.copy(postscript, file="figure1.eps", width=6, height=6, horizontal=F, onefile=F)
dev.off()

table(dat.subset$Is.the.microarray.data.publicly.available)
windows()
boxplot(Number.of.Citations.during.2004.2005 ~ Is.the.microarray.data.publicly.available,
        data = dat.subset,
        boxwex = 0.5, 
        names=c("Data Not Shared (n=43)", "Data Shared (n=27)"), 
        ylab = "Number of Citations in 2004-2005", outline=T, notch=F, log="y")
dev.copy(postscript, file="figure2.eps", width=6, height=6, horizontal=F, onefile=F)
dev.off()


