## setwd() to your folder where the file "Data_groom" is


# Loading packages --------------------------------------------------------
library("openxlsx") # version (‘4.2.5.2’)
library("DHARMa") # version (‘0.4.6’)
library("lme4") # version (‘1.1.35.3’)
library("car") # version (‘3.1.2’)
library(afex) # version (‘1.3.1’)
library(lmerTest) # version (‘3.1.3’)
library(glmmTMB) # version (‘1.1.9’)
library(sjPlot) # version (‘2.8.16’)
library(sjmisc) # version (‘2.8.10’)
library(sjlabelled) # version (‘1.2.0’)
library("performance") # version (‘0.12.0’)
library(plotrix) # version (‘3.8.4’)
library(ggplot2) # version (‘3.5.1’)
library('dplyr') # version (‘1.1.4’)
library(viridisLite) # version (‘0.4.2’)

# Loading data ------------------------------------------------------------
table = read.xlsx("Data_grooming.xlsx",colNames = TRUE) # Loading data file
table$Ind1[table$Ind1=="yvan"] = "yva"
table$Ind2[table$Ind2=="yvan"] = "yva"
## Transforming string / numerical data to factor ##
table$species_fac = as.factor(table$species)
table$temps_ECP = table$cou + table$torse + table$face 


rhesus_2021 = table[table$species== "Rhesus" & table$year==2021,]
tonk_2021 =table[table$species== "Tonkean" & table$year==2021,]
rhesus_2022 = table[table$species== "Rhesus" & table$year==2022,]
tonk_2023 = table[table$species== "Tonkean" & table$year==2023,]

# Calculate the z-score of the DSI, rank difference, age difference for each group
z.Rankdiff = rbind(scale(rhesus_2021$Rankdiff),scale(tonk_2021$Rankdiff),scale(rhesus_2022$Rankdiff),scale(tonk_2023$Rankdiff))
z.agediff = rbind(scale(rhesus_2021$agediff),scale(tonk_2021$agediff),scale(rhesus_2022$agediff),scale(tonk_2023$agediff))
z.Kin_degree = rbind(scale(rhesus_2021$kinship_degree),scale(tonk_2021$kinship_degree),scale(rhesus_2022$kinship_degree),scale(tonk_2023$kinship_degree))

#Create the column in the data table

table$z.agediff = as.vector(z.agediff)
table$z.Rankdiff = as.vector(z.Rankdiff)
table$male_pres_fac[table$male_pres==1]=TRUE
table$male_pres_fac[table$male_pres==0]=FALSE
table$z.Kin_degree = as.vector(z.Kin_degree)



table$Species = table$species_fac
table$Kin_degree = table$z.Kin_degree
table$Kin_degree_fac = table$kinship_degree
table$Rank_difference = table$z.Rankdiff
table$Age_difference = table$z.agediff
table$Male_presence = table$male_pres_fac
table$PS = table$prox_CP_EP / 210
table$reciprocity = table$reciprocite
## Creation of the table with the data used for the GLMM
my_data = table[, c('Species','Ind1','Ind2','Rank_difference','temps_toilettage','Male_presence',"Age_difference",'Kin_degree',"groomdos","groomface","temps_ECP","PS","reciprocity","groomflanc","groombd","Kin_degree_fac")]
my_data$temps_toilettage[my_data$temps_toilettage<30] = 0


# Proximity score analysis - Table 2 ------------------------------------------------

my_data$PS = (my_data$PS*(nrow(my_data) - 1) + 0.5)/nrow(my_data)
mdl_glm = glmmTMB(PS~1 +Species*(Kin_degree+Rank_difference)+Age_difference +Male_presence +(1|Ind1)+(1|Ind2) , data= my_data,family=beta_family)
summary(mdl_glm)

source("diagnostic_fcns.r")
overdisp.test(mdl_glm)
ranef.diagn.plot(mdl_glm)


source("glmmTMB_stability.r")
full.stab=glmmTMB.stab(model.res=mdl_glm, para=T,
                       data=my_data)

table(full.stab$detailed$converged)
m.stab.plot(full.stab$summary[, -1])

check_collinearity(mdl_glm)


tab_model(mdl_glm,transform = NULL,show.ci=F,show.se=T,show.re.var = F,  show.ngroups = F,show.r2=F,show.icc=F, string.se ="SE",show.obs = F,
          pred.labels = c("Intercept", "Species (Tonkean)", "Kin degree", "Rank difference",
                          "Age difference", "Presence of Male", 
                          "Tonkean : Kin degree","Tonkean : Rank difference"))

# Grooming time analysis - Table 3 / Figure 3-------------------------------------------------------
my_data$Involved_in_grooming[my_data$temps_toilettage<30] = 0 
my_data$Involved_in_grooming[my_data$temps_toilettage>30] = 1
my_data_without_zero = subset(my_data,my_data$temps_toilettage>30) #as we are looking for a proportion we deleted the dyad who groom less than 30 sec
my_data_without_zero$GT = my_data_without_zero$temps_toilettage

# analysis for 0-1 data
mdl_glm_binomial = glmmTMB( Involved_in_grooming~1 +Species*(Kin_degree+Rank_difference)+Age_difference +Male_presence+ (1|Ind1)+ (1|Ind2), data= my_data,family=binomial)
summary(mdl_glm_binomial)
simulateResiduals(mdl_glm_binomial, plot = T) # Residual (everything is ok)
check_collinearity(mdl_glm_binomial)
# analysis for data more than 0
mdl_glm = lmer( log(GT)~1 +Species*(Kin_degree+Rank_difference)+Age_difference +Male_presence+ (1|Ind1)+ (1|Ind2), data= my_data_without_zero)
summary(mdl_glm)
simulateResiduals(mdl_glm, plot = T) # Residual (everything is ok)
check_collinearity(mdl_glm)


tab_model(mdl_glm_binomial,mdl_glm,transform = NULL,show.ci=F,show.se=T,show.re.var = F,  show.ngroups = F,show.r2=F,show.icc=F, string.se ="SE",show.obs = F,string.est ="Estimates",
          pred.labels = c("Intercept", "Species (Tonkean)", "Kin degree", "Rank difference",
                          "Age difference", "Presence of Male", 
                          "Tonkean : Kin degree","Tonkean : Rank difference"),
          dv.labels = c("Involved in grooming", "log(GT)"))

### Graphical representation ###

### Rank difference 

ggplot(my_data, aes(as.factor(Involved_in_grooming),Rank_difference)) +
  geom_boxplot(width = 0.3)+
  scale_color_manual(values = c("grey", "black"))+
  theme_classic()+ theme(legend.position = "none")

### Kin degree*species

g <- my_data %>%
  group_by(Species,Kin_degree_fac,Involved_in_grooming) %>%
  summarize(cnt=n()) %>%
  mutate(mean = cnt / sum(cnt))%>%
  mutate(se = sqrt(mean*(1-mean)/cnt))

g = subset(g,g$Involved_in_grooming==1)

ggplot(data=g, aes(x=as.factor(Kin_degree_fac), y=mean, fill=Species)) + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_bar(stat='identity',position=position_dodge()) +
  geom_errorbar(
    stat='identity',
    width=0.2,
    aes(
      ymin=mean-se, 
      ymax = mean +se),position=position_dodge(.9)
  ) +
  scale_fill_manual(values = c("orange", "purple"))+
  theme_classic()

### Presence of male

# prepare a summary table with one row per experiment
df <- tibble(
  experiment=factor(c( rep('TRUE',length(my_data$Involved_in_grooming[my_data$Male_presence =='TRUE'])),
                       rep('FALSE', length(my_data$Involved_in_grooming[my_data$Male_presence =='FALSE'])) )),
  outcome=c( my_data$Involved_in_grooming[my_data$Male_presence =='TRUE'],my_data$Involved_in_grooming[my_data$Male_presence =='FALSE'])) %>%
  group_by(experiment) %>% 
  summarize(n=n(), p=mean(as.numeric(as.character(outcome)))) %>%
  mutate(se=sqrt(p*(1-p)/n))
colnames(df) <- c("Male_pres", "n","p","se")

# visualizing with barplot + errorbar
ggplot(data=df, aes(x=Male_pres, y=p, fill=Male_pres)) + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_bar(stat='identity') +
  geom_errorbar(
    stat='identity',
    width=0.5,
    aes(
      ymin=p-se, 
      ymax = p + se)
  ) +
  scale_fill_manual(values = c("grey47", "grey28"))+
  theme_classic()


## Kinship degree

ggplot(my_data_without_zero) +
  geom_boxplot(aes(as.factor(Kin_degree_fac),log(GT),colour = Species))+
  scale_color_manual(values = c("orange", "purple"))+
  theme_classic()+ theme(legend.position = "none")


# Face-to-back grooming analysis - Table 4 / Figure 4 ------------------------------------------
my_data_groom = subset(my_data,my_data$temps_toilettage>30) 
my_data_groom$Involved_in_fb_grooming[my_data_groom$groomdos==0] = 0 
my_data_groom$Involved_in_fb_grooming[my_data_groom$groomdos>0] = 1

mdl_glm_binomial = glmmTMB(Involved_in_fb_grooming~1 +Species*(Kin_degree+Rank_difference)+Age_difference +Male_presence + (1|Ind1)+ (1|Ind2), data= my_data_groom,family=binomial)

summary(mdl_glm_binomial)
simulateResiduals(mdl_glm_binomial, plot = T) # Residual (everything is ok)

check_collinearity(mdl_glm_binomial)

my_data_without_zero = subset(my_data_groom,my_data_groom$groomdos>0) #as we are looking for a proportion we deleted the dyad who groom less than 30 sec
my_data_without_zero$Pfb = my_data_without_zero$groomdos/my_data_without_zero$temps_toilettage

my_data_without_zero$Pfb = (my_data_without_zero$Pfb*(nrow(my_data_without_zero) - 1) + 0.5)/nrow(my_data_without_zero)


mdl_glm = glmmTMB(Pfb~1 +Species*(Kin_degree+Rank_difference)+Age_difference +Male_presence +(1|Ind1)+(1|Ind2) , data= my_data_without_zero,family=beta_family)
summary(mdl_glm)



source("diagnostic_fcns.r")
overdisp.test(mdl_glm)
ranef.diagn.plot(mdl_glm)
source("glmmTMB_stability.r")
full.stab=glmmTMB.stab(model.res=mdl_glm, para=T,
                       data=my_data_without_zero)

table(full.stab$detailed$converged)
m.stab.plot(full.stab$summary[, -1])

check_collinearity(mdl_glm)

tab_model(mdl_glm_binomial,mdl_glm,transform = NULL,show.ci=F,show.se=T,show.re.var = F,  show.ngroups = F,show.r2=F,show.icc=F, string.se ="SE",show.obs = F,string.est ="Estimates",
          pred.labels = c("Intercept", "Species (Tonkean)", "Kin degree", "Rank difference",
                          "Age difference", "Presence of Male", 
                          "Tonkean : Kin degree","Tonkean : Rank difference"),
          dv.labels = c("Involved in fb grooming", "Pfb"))
### Graphical representation ###

df <- tibble(
  experiment=factor(c( rep('Rhesus',length(my_data_groom$Involved_in_fb_grooming[my_data_groom$Species =='Rhesus'])),
                       rep('Tonkean', length(my_data_groom$Involved_in_fb_grooming[my_data_groom$Species =='Tonkean'])) )),
  outcome=c( my_data_groom$Involved_in_fb_grooming[my_data_groom$Species =='Rhesus'],my_data_groom$Involved_in_fb_grooming[my_data_groom$Species =='Tonkean'])) %>%
  group_by(experiment) %>% 
  summarize(n=n(), p=mean(as.numeric(as.character(outcome)))) %>%
  mutate(se=sqrt(p*(1-p)/n))
colnames(df) <- c("species", "n","p","se")



# visualizing with barplot + errorbar
ggplot(data=df, aes(x=species, y=p, fill=species)) + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_bar(stat='identity') +
  geom_errorbar(
    stat='identity',
    width=0.5,
    aes(
      ymin=p-se, 
      ymax = p + se)
  ) +
  scale_fill_manual(values = c("orange", "purple"))+
  theme_classic()

## Species


ggplot(my_data_without_zero,aes(Species,Pfb,colour = Species)) +
  geom_boxplot()+
  scale_color_manual(values = c("orange", "purple"))+
  geom_point()+
  geom_jitter(width=0.1)+
  theme_classic()+ theme(legend.position = "none")
# Kinship

ggplot(my_data_without_zero) +
  geom_boxplot(aes(as.factor(Kin_degree_fac),Pfb,colour = Species))+
  scale_color_manual(values = c("orange", "purple"))+
  theme_classic()+ theme(legend.position = "none")




# Face-to-face grooming analysis - Table 5 / Figure 5 ------------------------------------------
my_data_groom = subset(my_data,my_data$temps_toilettage>30) 

my_data_groom$Involved_in_ff_grooming[my_data_groom$groomface==0] = 0 
my_data_groom$Involved_in_ff_grooming[my_data_groom$groomface>0] = 1

mdl_glm_binomial = glmmTMB(Involved_in_ff_grooming~1 +Species*(Kin_degree+Rank_difference)+Age_difference +Male_presence+ (1|Ind1)+ (1|Ind2), data= my_data_groom,family=binomial)

summary(mdl_glm_binomial)
simulateResiduals(mdl_glm_binomial, plot = T) # Residual (everything is ok)

my_data_without_zero = subset(my_data_groom,my_data_groom$groomface>0) #as we are looking for a proportion we deleted the dyad who groom less than 30 sec
my_data_without_zero$Pff = my_data_without_zero$groomface/my_data_without_zero$temps_toilettage
my_data_without_zero$Pff = (my_data_without_zero$Pff*(nrow(my_data_without_zero) - 1) + 0.5)/nrow(my_data_without_zero)


mdl_glm = glmmTMB(Pff~1 +Species*(Kin_degree+Rank_difference)+Age_difference +Male_presence +(1|Ind1)+(1|Ind2) , data= my_data_without_zero,family=beta_family)
summary(mdl_glm)
source("diagnostic_fcns.r")
overdisp.test(mdl_glm)
ranef.diagn.plot(mdl_glm)


source("glmmTMB_stability.r")
full.stab=glmmTMB.stab(model.res=mdl_glm, para=T,
                       data=my_data_without_zero)

table(full.stab$detailed$converged)
m.stab.plot(full.stab$summary[, -1])
check_collinearity(mdl_glm)



tab_model(mdl_glm_binomial,mdl_glm,transform = NULL,show.ci=F,show.se=T,show.re.var = F,  show.ngroups = F,show.r2=F,show.icc=F, string.se ="SE",show.obs = F,string.est ="Estimates",
          pred.labels = c("Intercept", "Species (Tonkean)", "Kin degree", "Rank difference",
                          "Age difference", "Presence of Male", 
                          "Tonkean : Kin degree","Tonkean : Rank difference"),
          dv.labels = c("Involved in ff grooming", "Pff"))


### Graphical representation ###
### Species difference
df <- tibble(
  experiment=factor(c( rep(1,length(my_data_groom$Involved_in_ff_grooming[my_data_groom$Species =='Rhesus'])) , rep(2, length(my_data_groom$Involved_in_ff_grooming[my_data_groom$Species =='Tonkean'])) )),
  outcome=c( my_data_groom$Involved_in_ff_grooming[my_data_groom$Species =='Rhesus'],my_data_groom$Involved_in_ff_grooming[my_data_groom$Species =='Tonkean'])) %>%
  group_by(experiment) %>% 
  summarize(n=n(), p=mean(as.numeric(as.character(outcome)))) %>%
  mutate(se=sqrt(p*(1-p)/n))
colnames(df) <- c("species", "n","p","se")

# visualizing with barplot + errorbar
ggplot(data=df, aes(x=species, y=p, fill=species)) + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_bar(stat='identity') +
  geom_errorbar(
    stat='identity',
    width=0.5,
    aes(
      ymin=p-se, 
      ymax = p + se)
  ) +
  scale_fill_manual(values = c("orange", "purple"))+
  theme_classic()

### Kin degree * Species

g <- my_data_groom %>%
  group_by(Species,Kin_degree_fac,Involved_in_ff_grooming) %>%
  summarize(cnt=n()) %>%
  mutate(mean = cnt / sum(cnt))%>%
  mutate(se = sqrt(mean*(1-mean)/cnt))

g = subset(g,g$Involved_in_ff_grooming==1)

ggplot(data=g, aes(x=as.factor(Kin_degree_fac), y=mean, fill=Species)) + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_bar(stat='identity',position=position_dodge()) +
  geom_errorbar(
    stat='identity',
    width=0.2,
    aes(
      ymin=mean-se, 
      ymax = mean +se),position=position_dodge(.9)
  ) +
  scale_fill_manual(values = c("orange", "purple"))+
  theme_classic()

# Grooming of high ECP body parts - Table 6 / Figure 6 -----------------------------------------
my_data_groom = subset(my_data,my_data$temps_toilettage>30) 

my_data_groom$Involved_in_ECP_grooming[my_data_groom$temps_ECP==0] = 0 
my_data_groom$Involved_in_ECP_grooming[my_data_groom$temps_ECP>0] = 1

mdl_glm_binomial = glmmTMB(Involved_in_ECP_grooming~1 +Species*(Kin_degree+Rank_difference)+Age_difference +Male_presence+ (1|Ind1)+ (1|Ind2), data= my_data_groom,family=binomial)
summary(mdl_glm_binomial)
simulateResiduals(mdl_glm_binomial, plot = T) # Residual (everything is ok)


my_data_without_zero = subset(my_data_groom,my_data_groom$temps_ECP>0) #as we are looking for a proportion we deleted the dyad who groom less than 30 sec
my_data_without_zero$PEcp = my_data_without_zero$temps_ECP/my_data_without_zero$temps_toilettage
my_data_without_zero$PEcp = (my_data_without_zero$PEcp*(nrow(my_data_without_zero) - 1) + 0.5)/nrow(my_data_without_zero)


mdl_glm = glmmTMB(PEcp~1 +Species*(Kin_degree+Rank_difference)+Age_difference +Male_presence +(1|Ind1)+(1|Ind2) , data= my_data_without_zero,family=beta_family)
summary(mdl_glm)
source("diagnostic_fcns.r")
overdisp.test(mdl_glm)
ranef.diagn.plot(mdl_glm)

source("glmmTMB_stability.r")
full.stab=glmmTMB.stab(model.res=mdl_glm, para=T,
                       data=my_data_without_zero)

table(full.stab$detailed$converged)
m.stab.plot(full.stab$summary[, -1])
xx=lme4::lmer(PEcp~1 +Species*(Kin_degree+Rank_difference)+Age_difference +Male_presence +(1|Ind1)+(1|Ind2) ,
              data=my_data_without_zero)
round(vif(xx), 3)


tab_model(mdl_glm_binomial,mdl_glm,transform = NULL,show.ci=F,show.se=T,show.re.var = F,  show.ngroups = F,show.r2=F,show.icc=F, string.se ="SE",show.obs = F,string.est ="Estimates",
          pred.labels = c("Intercept", "Species (Tonkean)", "Kin degree", "Rank difference",
                          "Age difference", "Presence of Male", 
                          "Tonkean : Kin degree","Tonkean : Rank difference"),
          dv.labels = c("Involved in ECP grooming", "Pecp"))
### Graphical representation ###
# Species


df <- tibble(
  experiment=factor(c( rep(1,length(my_data_groom$Involved_in_ECP_grooming[my_data_groom$Species =='Rhesus'])) , rep(2, length(my_data_groom$Involved_in_ECP_grooming[my_data_groom$Species =='Tonkean'])) )),
  outcome=c( my_data_groom$Involved_in_ECP_grooming[my_data_groom$Species =='Rhesus'],my_data_groom$Involved_in_ECP_grooming[my_data_groom$Species =='Tonkean'])) %>%
  group_by(experiment) %>% 
  summarize(n=n(), p=mean(as.numeric(as.character(outcome)))) %>%
  mutate(se=sqrt(p*(1-p)/n))
colnames(df) <- c("species", "n","p","se")

# visualizing with barplot + errorbar
ggplot(data=df, aes(x=species, y=p, fill=species)) + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_bar(stat='identity') +
  geom_errorbar(
    stat='identity',
    width=0.5,
    aes(
      ymin=p-se, 
      ymax = p + se)
  ) +
  scale_fill_manual(values = c("orange", "purple"))+
  theme_classic()

# Male presence

df <- tibble(
  experiment=factor(c( rep('TRUE',length(my_data_groom$Involved_in_ECP_grooming[my_data_groom$Male_presence =='TRUE'])),
                       rep('FALSE', length(my_data_groom$Involved_in_ECP_grooming[my_data_groom$Male_presence =='FALSE'])) )),
  outcome=c( my_data_groom$Involved_in_ECP_grooming[my_data_groom$Male_presence =='TRUE'],my_data_groom$Involved_in_ECP_grooming[my_data_groom$Male_presence =='FALSE'])) %>%
  group_by(experiment) %>% 
  summarize(n=n(), p=mean(as.numeric(as.character(outcome)))) %>%
  mutate(se=sqrt(p*(1-p)/n))
colnames(df) <- c("Male_pres", "n","p","se")



ggplot(data=df, aes(x=Male_pres, y=p, fill=Male_pres)) + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_bar(stat='identity') +
  geom_errorbar(
    stat='identity',
    width=0.5,
    aes(
      ymin=p-se, 
      ymax = p + se)
  ) +
  scale_fill_manual(values = c("grey47", "grey28"))+
  theme_classic()
### Kin degree * Species

g <- my_data_groom %>%
  group_by(Species,Kin_degree_fac,Involved_in_ECP_grooming) %>%
  summarize(cnt=n()) %>%
  mutate(mean = cnt / sum(cnt))%>%
  mutate(se = sqrt(mean*(1-mean)/cnt))

g = subset(g,g$Involved_in_ECP_grooming==1)

ggplot(data=g, aes(x=as.factor(Kin_degree_fac), y=mean, fill=Species)) + 
  coord_cartesian(ylim=c(0,1)) + 
  geom_bar(stat='identity',position=position_dodge()) +
  geom_errorbar(
    stat='identity',
    width=0.2,
    aes(
      ymin=mean-se, 
      ymax = mean +se),position=position_dodge(.9)
  ) +
  scale_fill_manual(values = c("orange", "purple"))+
  theme_classic()

# Age difference 

pt.cols=plasma(n=2, begin=0.8, end=0.2, alpha=0.5)
l.cols=plasma(n=2, begin=0.8, end=0.2)
coefs=coefs=fixef(mdl_glm)$cond

par(mar=c(3, 3, 0.2, 0.2), mgp=c(1.7, 0.3, 0), tcl=-0.15,
    las=1)
plot(x=my_data_without_zero$Age_difference, y=my_data_without_zero$PEcp,
     pch=19, col=pt.cols[as.numeric(my_data_without_zero$Species)],
     xlab="Normalized age difference", ylab="Proportion of grooming of high ECP bodypats",
     ylim=c(0, 1), xaxt="n")
xlab=seq(from=-2, to=10, by=0.5)
#map to z-space:
xat=(xlab-mean(my_data_without_zero$Age_difference))/sd(my_data_without_zero$Age_difference)
#add axis:
axis(side=1, at=xat, labels=xlab)

xvals=seq(from=min(my_data_without_zero$Age_difference), to=max(my_data_without_zero$Age_difference),
          length.out=100)
#correct condition:
yvals.corr=coefs["(Intercept)"]+coefs["Age_difference"]*xvals
#incorrect condition:
yvals.incorr=coefs["(Intercept)"] + coefs["SpeciesTonkean"] +
  (coefs["Age_difference"] )*xvals

#transform to response space:
yvals.corr=exp(yvals.corr)/(1+exp(yvals.corr))
yvals.incorr=exp(yvals.incorr)/(1+exp(yvals.incorr))
lines(x=xvals, y=yvals.incorr, col=l.cols[2], lwd=3, lty=2)
lines(x=xvals, y=yvals.corr, col=l.cols[1], lwd=3, lty=2)

# Supplementary - Reciprocity ---------------------------------------------


my_data_groom = subset(my_data,my_data$temps_toilettage>30) 
my_data_groom$zi_rec[my_data_groom$reciprocity==0] = 0 
my_data_groom$zi_rec[my_data_groom$reciprocity>0] = 1

mdl_glm_binomial = glmmTMB(zi_rec~1 +Species*(Kin_degree+Rank_difference)+Age_difference +Male_presence + (1|Ind1)+ (1|Ind2), data= my_data_groom,family=binomial)

summary(mdl_glm_binomial)
simulateResiduals(mdl_glm_binomial, plot = T) # Residual (everything is ok)

check_collinearity(mdl_glm_binomial)

my_data_without_zero = subset(my_data_groom,my_data_groom$reciprocity>0) #as we are looking for a proportion we deleted the dyad who groom less than 30 sec

mdl_glm = glmmTMB(reciprocity~1 +Species*(Kin_degree+Rank_difference)+Age_difference +Male_presence +(1|Ind1)+(1|Ind2) , data= my_data_without_zero,family=beta_family)
summary(mdl_glm)



source("diagnostic_fcns.r")
overdisp.test(mdl_glm)
ranef.diagn.plot(mdl_glm)
source("glmmTMB_stability.r")
full.stab=glmmTMB.stab(model.res=mdl_glm, para=T,
                       data=my_data_without_zero)

table(full.stab$detailed$converged)
m.stab.plot(full.stab$summary[, -1])

check_collinearity(mdl_glm)

tab_model(mdl_glm_binomial,mdl_glm,transform = NULL,show.ci=F,show.se=T,show.re.var = F,  show.ngroups = F,show.r2=F,show.icc=F, string.se ="SE",show.obs = F,string.est ="Estimates",
          pred.labels = c("Intercept", "Species (Tonkean)", "Kin degree", "Rank difference",
                          "Age difference", "Presence of Male", 
                          "Tonkean : Kin degree","Tonkean : Rank difference"),
          dv.labels = c("Involved in reciprocal grooming", "Reciprocity"))



