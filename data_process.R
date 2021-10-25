library("rstan")
library("bayesplot")
library("zoo")
library("dplyr")
library("reshape2")
library("scales")
library("ggridges")
library("ggthemes")
library("ggcor")
library("ggplot2")
library("gridExtra")
library("ggalt")
library("aomisc")
library("meta")
library("RColorBrewer")
library("grid")
setwd("E:/COVID_Vaccination/NPI&Vaccination_code")#set up your setwd here
Sys.setlocale("LC_TIME","English")
source('main.R')
source('figurecode.R')

RT<-read.csv("Dataset/owid-covid-data.csv")# Daily reproduction rate data (Rt)
Variants<-read.csv("Dataset/covid-variants.csv")# Weekly variant structure data
NPI<-"Dataset/OxCGRT_latest_0923.csv" # Daily non-pharmaceutical interventions data (NPI)
POP<-"Dataset/pop_data_age_by_country.csv" # Total population number of each country
index<-'Dataset/index.csv'# Control variables, contained population density, age structre and health index of each country
Env<-'Dataset/20200120-20210923humidity&airtem.csv'# Control variables, contained temperature and humidity. Only air temperature were considered in our study
Vaccination<-'Dataset/vaccinations.csv' # Daily vaccination data 
Vaccman<-'Dataset/vaccinations-by-manufacturer.csv'# Weekly vaccination manufacturer structure

school_holiday<-"Dataset/day_public_and_school_holidays_2010_2019.csv"
variant_tran<-"Dataset/variant_transmission.csv" #transmission parameters of each variant
vaccine_effect<-"Dataset/vaccine_effectiveness_V1.csv" #V1:the default setting 
                                                       #S1:the absent effects of vaccines were the same as that against Alpha
                                                       #S2:the absent effects of vaccines were the same as that against Delta

####################################################################################################
########### part I: Effect of NPIs and vaccination over time on country level ###########

path<-"SIresult_forAlpha" #set up your output path here
for(i in seq(1,30)){
  if (dir.exists(path)==F){dir.create(path)
  }else{print("This path has been exists")}
  dataprocess(RT,Variants,Vaccination,NPI,POP,index,Env,
              vaccine_effect,variant_tran,google_mobility,i,path)
  gc()
}
figpath<-"0924picture"# set up the output path of the figures 
country_result_plot(path,figpath)

#########################################################################################

############# part II: regional result based on Meta-analysis ###########################

meta_analysis(path)

###################################################################################################

###################### part III: Impact of vaccination on the effectiveness of NPIs ###############

library("rcompanion")
library("tidyverse")
library("investr")
library("modelr")
library("grid")
library("ggpubr")

fig3plot(figpath)

##################################################################################################################
#################################part IV: Relaxation of NPIs amid vaccination################################

library("rstatix")
risk<-read.csv("dataset/riskindex_timeseries_latest.csv",stringsAsFactors = F)# an openness risk calculated by Hale et al.
r<-fig4plot(risk,figpath)
r%>%select(openness_risk,Open)%>%cor_test()
cor.test(r$openness_risk,r$Open, type="pearson") 

#############################################################################################################
####RSS and Rhat#########

McmcParplot(path)

###########################################
#############Correlation analysis##########

Dataset<-read.csv("dataset/Europe_dataset_0923V1.csv",row.names = F)
N<-Dataset[,c("R0","Rt","StringencyIndex","Fully_vaccinated_effect","Tem","Hum")]
colnames(N)<-c("R0","R0,t","Stringency index","Efficient vaccinated rate","Temperature","Humidity")
cor<-COR(N)
ggsave(paste0(figpath,"/correlation_analysis.pdf"),cor,width=120,height=140,units="mm",device = cairo_pdf)

###############################################################################################################
###############Figure 1 plot##############

Dataset<-read.csv("dataset/Europe_dataset_0923V1.csv",row.names = F)
fig1a1plot(Dataset,figpath)
fig1a2plot(Dataset,figpath)
fig1bplot(Dataset,figpath)
fig1cplot(Dataset,figpath)
fig1dplot(Dataset,figpath)

###############Figure 2 plot##############

library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("dplyr")
library("mapproj")
library("fiftystater")

fig2a1plot(figpath)
fig2a2plot(Dataset,figpath)

#################################################################################################################


