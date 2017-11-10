library(PheWAS)
## PheWAS with a phecode based propensity score
## Create original data objects
### Demographic data
#### Input should contain 'ID', 'AGE', 'RACE', 'SEX'
#### Our race variable is binary for EHR recorded 'White' versus all others.
#### Our sex variable is binary for EHR recorded 'Female'. Unknown individuals were excluded.
demo=demo.raw %>% select(ID, AGE, RACE, SEX)
### Phecodes
#### Input data should be of the form 'ID', 'ICD9', 'COUNT' - note that count is not really used here
#### Baseline phecodes were generated from all ICD-9-CM codes where age was 0 years

phewas_codes.expanded=rbind(phewas_codes,data.frame(ID=demo$ID,PHEWAS_CODE="DUMMY",COUNT=10))
phewas_codes_baseline.expanded=rbind(phewas_codes_baseline,data.frame(ID=demo$ID,PHEWAS_CODE="DUMMY",COUNT=10))
pheno=createPhewasTable(phewas_codes.expanded,translate = F, min.code.count=1, #id.gender = demo.raw %>% select(ID,GENDER_EPIC),
                        add.exclusions = F) %>% select(-DUMMY)
pheno.baseline=createPhewasTable(phewas_codes_baseline.expanded,translate = F,min.code.count = 1,add.exclusions = F) %>% select(-DUMMY)
exposure = exposure.raw %>% mutate(EXPOSURE=EXPOSURE=='t') %>% filter(!duplicated(ID))
d=inner_join(inner_join(exposure,demo),pheno)
d.baseline=inner_join(inner_join(exposure,demo),pheno.baseline)

d=inner_join(id_map,d) %>% select(-ID)
d.baseline=inner_join(id_map,d.baseline) %>% select(-ID)

#Process the data

d.baseline.org <- d.baseline %>% select(-GENDER_EPIC,-SD_RACE_ETH)

d.org <- d %>% select(-GENDER_EPIC,-SD_RACE_ETH)

#### Outcome phenotype data preprocessing
#Select only those phenotypes with more than one case
all.phenotypes <- colSums(d.org[,-1:-5],na.rm = T)
phenotype.final.list <- names(all.phenotypes[all.phenotypes>1]) # exclude all false and only one case

#Normalize age
dd.O.all <- subset(d.org, select=c('GRID','EXPOSURE','AGE','SD_RACE_W','GENDER_EPIC_F', phenotype.final.list))
dd.O.all$AGE <- (dd.O.all$AGE - mean(dd.O.all$AGE))/sd(dd.O.all$AGE)

#### covariate (baseline) data preprocessing
all.covariates <- colnames(d.baseline.org)[-1:-5]
#Select only 3digit phecodes
all.covariates <- all.covariates[nchar(all.covariates)==3]
#Remove those with no cases
all.covariates.sums = colSums(d.baseline.org[,all.covariates],na.rm = T)
covariates.baseline.all <- names(all.covariates.sums[all.covariates.sums>0])

dd.all <- d.baseline.org[, c('GRID','EXPOSURE','AGE','SD_RACE_W','GENDER_EPIC_F', covariates.baseline.all)]
dd.all$AGE <- (dd.all$AGE - mean(dd.all$AGE))/sd(dd.all$AGE)

#### predict PS and generate the final outcome data for PheWAS
set.seed(10)
alpha.para <- 0.1
data.x <- as.matrix(dd.all[,-c(1:2)])
glmnet.fit <- cv.glmnet(x=data.x, y=dd.all$EXPOSURE, family="binomial", standardize=TRUE, alpha=alpha.para)
ps <- predict(glmnet.fit, data.x, s='lambda.min')[,"1"]
data.ps <- data.frame(GRID=dd.all$GRID, PS=ps)
dd.all.ps <- inner_join(data.ps, dd.O.all)


####Setup PheWAS analysis
dd.all.ps[1:10,1:10]

res=phewas_ext(phenotypes=phenotype.final.list,predictors="EXPOSURE",covariates="PS",
               data=dd.all.ps,additive.genotypes=F,method="logistf",min.records=0, cores = 40)

amp.res=res %>% mutate(predictor="Ampicillin")

save(amp.res,file="~/Dropbox/research/kids/matching/amp.mcc1.noex.RData")

res.demo=phewas_ext(phenotypes=phenotype.final.list,predictors="EXPOSURE",covariates=c('AGE','SD_RACE_W','GENDER_EPIC_F'),
                    data=dd.all.ps,additive.genotypes=F,method="logistf",min.records=0, cores = 40)

res.demo=res.demo %>% mutate(predictor="Ampicillin")

save(res.demo,file="~/Dropbox/research/kids/matching/amp.mcc1.noex.demographics.RData")


#### Plot Results

load("~/Dropbox/research/kids/matching/amp.mcc1.noex.RData")

#res.pml.amp=res.pml.amp %>% filter(n.cases>=20)
amp= amp.res %>%  addPhecodeInfo() %>% filter(group!='congenital anomalies',!grepl("congenital",description,ignore.case=T)) %>% arrange(p)

a.plot=phewasManhattan(amp %>% transmute(phenotype=phecode,OR,p,n_cases),
                       annotate.level=.05,title="PheWAS on ampamicin Exposure, PS Adjusted",y.axis.interval=1,OR.direction=T,
                       significant.line=NA, max.y=5) +
  aes(size=cut(n_cases,breaks=c(0,20,50,100,Inf),right=FALSE)) + scale_size_discrete(range=c(2,6),guide="legend") + guides(size="legend")
ggsave(a.plot,filename = "~/Dropbox/research/kids/matching/amp.mcc1.noex.png",width = 8, height=4.5,scale = 2.2)
amp %>% head()

amp %>% write.csv(file="~/Dropbox/research/kids/matching/amp.mcc1.noex.csv",row.names=F)


load("~/Dropbox/research/kids/matching/amp.mcc1.noex.demographics.RData")

#res.pml.amp=res.pml.amp %>% filter(n.cases>=20)
amp.demo = res.demo %>%  addPhecodeInfo() %>% filter(group!='congenital anomalies',!grepl("congenital",description,ignore.case=T)) %>% arrange(p)

a.d.plot=phewasManhattan(amp.demo %>% transmute(phenotype=phecode,OR,p,n_cases),
                         annotate.level=1e-14,title="PheWAS on ampamicin Exposure, demographics Adjusted",OR.direction=T) +
  aes(size=cut(n_cases,breaks=c(0,20,50,100,Inf),right=FALSE)) + scale_size_discrete(range=c(2,6),guide="legend") + guides(size="legend")
ggsave(a.d.plot,filename = "~/Dropbox/research/kids/matching/amp.mcc1.noex.demo.png",width = 8, height=4.5,scale = 2.2)
amp.demo %>% head()

amp.demo %>% write.csv(file="~/Dropbox/research/kids/matching/amp.mcc1.noex.demo.csv",row.names=F)
