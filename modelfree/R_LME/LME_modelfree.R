library(readr)
library(lmerTest)
library(sjPlot)
df<- read.csv("D:/Trust_Game/Data_summary/mean_tiles_covered_TD.csv")
df<- read.csv("D:/Trust_Game/Data_summary/mean_tiles_covered_TD_matched.csv")

df[,'sage']<-scale(df[,'age'], center = TRUE, scale = TRUE)
df[,'srs']<-scale(df[,'srs'], center = TRUE, scale = TRUE)
df[,'subject']<-factor(df[,'subject'])
df[,'gender']<-factor(df[,'gender'])



###infosample###
# gender
modela <- lmer(sampled ~ gender + (1|subject), data=df)
tab_model(modela)
summary(modela)


###############################################################################################
#TD uncertainty
model1 <- lmer(sampled ~ uncertainty + (1|subject), data=df)
tab_model(model1)
summary(model1)

#TD uncertainty age
model2 <- lmer(sampled ~ uncertainty*poly(sage,degree = 2,raw = TRUE)+gender + (1|subject), data=df)
tab_model(model2)
summary(model2)

model3 <- lmer(sampled ~ uncertainty*srs + (1|subject), data=df)
tab_model(model3)
summary(model3)


model4 <- lmer(sampled ~ uncertainty*poly(sage,degree = 2,raw = TRUE)*srs + (1|subject), data=df)
tab_model(model4)
summary(model4)
############################################################################################
#TD reciprocation
model11 <- lmer(investprob ~ reciprocation + (1|subject), data=df)
tab_model(model11)
summary(model11)

#TD reciprocation age
model12 <- lmer(investprob ~ reciprocation*poly(sage,degree = 2,raw = TRUE) + (1|subject), data=df)
tab_model(model12)
summary(model12)

#TD reciprocation age 
model3 <- lmer(investprob ~ reciprocation*srs + (1|subject), data=df)
tab_model(model3)
summary(model3)

model14 <- lmer(investprob ~ reciprocation*sage*srs + (1|subject), data=df)
tab_model(model14)
summary(model14)
#################################################################################
#################################################################################
df<- read.csv("D:/Trust_Game/Data_summary/mean_tiles_covered_all_matched.csv")

df[,'sage']<-scale(df[,'age'], center = TRUE, scale = TRUE)
df[,'srs']<-scale(df[,'srs'], center = TRUE, scale = TRUE)
df[,'develop']<-factor(df[,'develop'], levels = c('TD','ASD'))
df[,'subject']<-factor(df[,'subject'])
df[,'gender']<-factor(df[,'gender'])


#TD uncertainty
model1 <- lmer(sampled ~ uncertainty*develop + (1|subject), data=df)
tab_model(model1)
summary(model1)

model2 <- lmer(sampled ~ uncertainty*poly(sage,degree = 2,raw = TRUE)*develop*srs + (1|subject), data=df)
tab_model(model2)
summary(model2)

model2 <- lmer(sampled ~ uncertainty*sage*develop*srs + (1|subject), data=df)
tab_model(model2)
summary(model2)
############################################################################################
#TD reciprocation
model11 <- lmer(investprob ~ reciprocation*develop +(1|subject), data=df)
tab_model(model11)
summary(model11)

#TD reciprocation age srs/aq
model13a <- lmer(investprob ~ reciprocation*sage*srs*develop*gender + (1|subject), data=df)
tab_model(model13a)
summary(model13a)
#srs
model13b <- lmer(investprob ~ reciprocation*sage*srs*develop + (1|subject), data=df)
tab_model(model13b)
summary(model13b)


#################################################################################
#################################################################################
#################################################################################
#################################################################################
library(tidyverse)
library(extrafont)
############################################################################################
# fig2 B D
df <- read.csv("D:/Online_Trust_Game/Data_files/All_descriptive_male.csv")
df <- read.csv("D:/Online_Trust_Game/Data_files/All_descriptive_male_young.csv")
df <- read.csv("D:/Online_Trust_Game/Data_files/All_descriptive_male_old.csv")

df <- read.csv("D:/Online_Trust_Game/Data_files/All_descriptive_matched.csv")
df <- read.csv("D:/Online_Trust_Game/Data_files/ASD_descriptive_matched.csv")
df <- read.csv("D:/Online_Trust_Game/Data_files/TD_descriptive_matched.csv")
#age invest
ggplot(df,aes(x=age,y=invest,color=Reciprocation))+
  geom_point(size=5)+
  geom_smooth(method="lm",aes(fill=Reciprocation),size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
  theme(axis.line.x = element_line(color = 'black',size = 1),
        axis.line.y = element_line(color = 'black',size = 1))+
  
  labs(fill="Reciprocation", x="Age(month)", y="Proportion of invest")+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand = c(0, 0), limits=c(-2,2))+
  coord_cartesian(ylim=c(-0.05,1.05),xlim=c(70,190))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))

scale_fill_manual(values=c("#999999", "#E69F00")) + 
  scale_color_manual(values=c("#999999", "#E69F00")) +
  #srs invest
  ggplot(df,aes(x=SRS,y=invest,color=Reciprocation))+
  geom_point(size=5)+
  geom_smooth(method="lm",aes(fill=Reciprocation),size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
  theme(axis.line.x = element_line(color = 'black',size = 1),
        axis.line.y = element_line(color = 'black',size = 1))+
  
  labs(fill="Reciprocation", x="SRS Raw", y="Proportion of invest")+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand = c(0, 0), limits=c(-2,2))+
  coord_cartesian(ylim=c(-0.05,1.05),xlim=c(0,180))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))

scale_fill_manual(values=c("#999999", "#E69F00")) + 
  scale_color_manual(values=c("#999999", "#E69F00")) +
############################################################################################
# fig3 C
df <- read.csv("D:/Online_Trust_Game/Data_files/ASD_descriptive_matched_young.csv")
df <- read.csv("D:/Online_Trust_Game/Data_files/TD_descriptive_matched_young.csv")  
df <- read.csv("D:/Online_Trust_Game/Data_files/ASD_descriptive_matched_middle.csv")
df <- read.csv("D:/Online_Trust_Game/Data_files/TD_descriptive_matched_middle.csv") 
df <- read.csv("D:/Online_Trust_Game/Data_files/ASD_descriptive_matched_old.csv")
df <- read.csv("D:/Online_Trust_Game/Data_files/TD_descriptive_matched_old.csv") 
#srs sample
ggplot(df2,aes(x=SRS,y=tiles,color=Uncertainty))+
  geom_point(size=5)+
  geom_smooth(method="lm",aes(fill=Uncertainty),size = 3)+
  theme_classic()+
  theme(text=element_text(size=30, family = "Arial",face="bold"))+
  theme(axis.line.x = element_line(color = 'black',size = 1),
        axis.line.y = element_line(color = 'black',size = 1))+
  labs(fill="Uncertainty", x="SRS", y="Num of tiles sampled")+
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand = c(0, 0), limits=c(-10,40))+
  coord_cartesian(ylim=c(-0.35,25.9),xlim=c(0,180))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+