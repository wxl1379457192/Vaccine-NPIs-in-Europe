Vaccination_complementation<-function(Vaccination){
  V1<-read.csv(Vaccination,stringsAsFactors = F)
  V1<-do.call(rbind,lapply(split(V1,V1$iso_code),FUN=function(x){
    x$date<-as.Date(x$date)
    startdate <- min(x$date)
    enddate <-max(x$date)
    date2 <- as.data.frame(seq.Date(from = startdate,to=enddate,by="days"))
    colnames(date2)<-"Date"
    date2$location<-unique(x$location)
    date2$iso_code<-unique(x$iso_code)
    date2$total_vaccinations<-0
    date2$people_vaccinated<-0
    date2$people_fully_vaccinated<-0
    date2$daily_vaccinations_raw<-0
    date2$daily_vaccinations<-0
    date2$total_vaccinations_per_hundred<-0
    date2$people_vaccinated_per_hundred<-0
    date2$people_fully_vaccinated_per_hundred<-0
    date2$daily_vaccinations_per_million<-0
    for(i in 4:ncol(x)){
      for(j in 1:nrow(x)){
        d<-x$date[j]
        date2[which(date2$Date==d),i]<-x[j,i]
      }
    }
    return(date2)
  }))
  #write.csv(V1,file = 'Dataset/vaccinations.csv',row.names = F)
  return(V1)
}
nafill<-function(x){
  if(length(which(is.na(x)))>0){
    list<-which(is.na(x))
    for (i in 1:length(list)){
      ii<-list[i]%%7
      l<-seq(ii,length(x),7)
      data<-zoo(x[l])
      x[l]<-na.fill(data,"extend")
    }
    return(x)
  }else(return(x))
}
R0calculate<-function(RT,Variants,R0_base){
  RT<-RT[,c("iso_code","continent","location","date","reproduction_rate")]
  RT$date<-as.Date(RT$date)
  #vartran<-read.csv("Dataset/variant_tranmission.csv")
  Variants[is.na(Variants)]<-0
  Variants$date<-as.Date(Variants$date)
  ROlist<-lapply(split(RT,RT$iso_code),FUN = function(v){
    V<-subset(Variants,Variants$location==unique(v$location))
    r0<-R0_base$R0[which(R0_base$location==unique(v$location))]
    if (nrow(V)>0&length(r0>0)){
      v$R0<-r0
      v$R0_low<-r0
      v$R0_upper<-r0
      datelist<-split(V,V$date)
      for (i in datelist){
        if(unique(i$date)-14>=min(v$date)){
          alphaN<- as.numeric(i$perc_sequences[which(i$variant=="Alpha")])/100
          betaN<-as.numeric(i$perc_sequences[which(i$variant=="Beta")])/100
          gammaN<-as.numeric(i$perc_sequences[which(i$variant=="Gamma")])/100
          deltaN<-as.numeric(i$perc_sequences[which(i$variant=="Delta")])/100
          etaN<-as.numeric(i$perc_sequences[which(i$variant=="Eta")])/100
          kappaN<-as.numeric(i$perc_sequences[which(i$variant=="Kappa")])/100
          baseN<-1-alphaN-betaN-gammaN-deltaN-etaN-kappaN
          v$R0[(which(v$date==unique(i$date))-14):which(v$date==unique(i$date))]<-
            r0*1.29*alphaN+r0*1.25*betaN+r0*1.38*gammaN+
            r0*1.97*deltaN+r0*1.29*etaN+r0*1.48*kappaN+r0*baseN
          v$R0_low[(which(v$date==unique(i$date))-14):which(v$date==unique(i$date))]<-
            r0*1.24*alphaN+r0*1.20*betaN+r0*1.29*gammaN+
            r0*1.76*deltaN+r0*1.23*etaN+r0*1.28*kappaN+r0*baseN
          v$R0_upper[(which(v$date==unique(i$date))-14):which(v$date==unique(i$date))]<-
            r0*1.33*alphaN+r0*1.30*betaN+r0*1.48*gammaN+
            r0*2.17*deltaN+r0*1.35*etaN+r0*1.69*kappaN+r0*baseN
        }
      }
      m <- v[complete.cases(v),]
      print(paste(unique(m$location),"has been processed"))
      return(m)
    }
    return(NULL)
  })
  R0data<-do.call(rbind,ROlist)
  return(R0data)
}
approx_sd<-function(x1,x2){
  (x2-x1)/(qnorm(0.95)-qnorm(0.05))
}
#vartran$sd<-approx_sd(vartran$IQR.,vartran$IQR..1)
variant_add<-function(RT,Variants){
  #combined the variant structure data and reproduction data
  RT<-RT[,c("iso_code","continent","location","date","reproduction_rate","new_cases","new_cases_smoothed")]
  RT$date<-as.Date(RT$date)
  Variants[is.na(Variants)]<-0
  Variants$date<-as.Date(Variants$date)
  colnames(RT)<-c("iso_code","continent","location","Date","Rt","Total_cases","New_cases_smoothed")
  list<-lapply(split(RT,RT$location),FUN = function(v){
    V<-subset(Variants,Variants$location==unique(v$location))
    if (nrow(V)>0){
      v$alphaN<-as.numeric(0)
      v$betaN<-as.numeric(0)
      v$gammaN<-as.numeric(0)
      v$etaN<-as.numeric(0)
      v$kappaN<-as.numeric(0)
      v$deltaN<-as.numeric(0)
      v$baseN<-as.numeric(1)
      datelist<-split(V,V$date)
      for (i in datelist){
        if(unique(i$date)-14>=min(v$Date)){
          v$alphaN[(which(v$Date==unique(i$date))-13):which(v$Date==unique(i$date))]<-
            as.numeric(i$perc_sequences[which(i$variant=="Alpha")])/100
          v$betaN[(which(v$Date==unique(i$date))-13):which(v$Date==unique(i$date))]<-
            as.numeric(i$perc_sequences[which(i$variant=="Beta")])/100
          v$gammaN[(which(v$Date==unique(i$date))-13):which(v$Date==unique(i$date))]<-
            as.numeric(i$perc_sequences[which(i$variant=="Gamma")])/100
          v$deltaN[(which(v$Date==unique(i$date))-13):which(v$Date==unique(i$date))]<-
            as.numeric(i$perc_sequences[which(i$variant=="Delta")])/100
          v$etaN[(which(v$Date==unique(i$date))-13):which(v$Date==unique(i$date))]<-
            as.numeric(i$perc_sequences[which(i$variant=="Eta")])/100
          v$kappaN[(which(v$Date==unique(i$date))-13):which(v$Date==unique(i$date))]<-
            as.numeric(i$perc_sequences[which(i$variant=="Kappa")])/100
          a<-1-as.numeric(i$perc_sequences[which(i$variant=="Alpha")])/100-
            as.numeric(i$perc_sequences[which(i$variant=="Beta")])/100-
            as.numeric(i$perc_sequences[which(i$variant=="Gamma")])/100-
            as.numeric(i$perc_sequences[which(i$variant=="Delta")])/100-
            as.numeric(i$perc_sequences[which(i$variant=="Eta")])/100-
            as.numeric(i$perc_sequences[which(i$variant=="Kappa")])/100
          ifelse(a>=0,v$baseN[(which(v$Date==unique(i$date))-13):which(v$Date==unique(i$date))]<-a,0)
        }
        if(max(v$Date)>max(V$date)){
          v$deltaN[which(v$Date>max(V$date))]<-v$deltaN[which(v$Date==max(V$date))]
          v$alphaN[which(v$Date>max(V$date))]<-v$alphaN[which(v$Date==max(V$date))]
          v$betaN[which(v$Date>max(V$date))]<-v$betaN[which(v$Date==max(V$date))]
          v$gammaN[which(v$Date>max(V$date))]<-v$gammaN[which(v$Date==max(V$date))]
          v$etaN[which(v$Date>max(V$date))]<-v$etaN[which(v$Date==max(V$date))]
          v$kappaN[which(v$Date>max(V$date))]<-v$kappaN[which(v$Date==max(V$date))]
          v$baseN[which(v$Date>max(V$date))]<-v$baseN[which(v$Date>max(V$date))]-
            v$alphaN[which(v$Date==max(V$date))]-v$betaN[which(v$Date==max(V$date))]-
            v$gammaN[which(v$Date==max(V$date))]-v$etaN[which(v$Date==max(V$date))]-
            v$kappaN[which(v$Date==max(V$date))]-v$deltaN[which(v$Date==max(V$date))]
        }
      }
      m <- v[complete.cases(v),]
      print(paste(unique(m$location),"has been processed"))
      return(m)
    }else{return(NULL)}
  })
  data<-do.call(rbind,list)
  return(data)
}
vaccine_ratio_cal<-function(Vman,V1){
  Vman$date<-as.Date(Vman$date)
  V1$Date<-as.Date(V1$Date)
  V1<-V1[,1:4]
  V1$total_vaccinations[is.na(V1$total_vaccinations)]<-0
  colnames(Vman)<-c("location","Date","vaccine","total_vaccinations")
  list<-lapply(split(Vman,Vman$location),FUN = function(v){
    vac<-subset(V1,V1$location==unique(v$location))
    vac$total_vaccinations<-cummax(vac$total_vaccinations)
    if (nrow(vac)>0){
      vac$'Johnson&Johnson'<-as.numeric(0)
      vac$'Moderna'<-as.numeric(0)
      vac$'Oxford/AstraZeneca'<-as.numeric(0)
      vac$'Pfizer/BioNTech'<-as.numeric(0)
      vac$'Sinovac'<-as.numeric(0)
      vac$'CanSino'<-as.numeric(0)
      vac$'Sputnik V'<-as.numeric(0)
      vac$'Sinopharm/Beijing'<-as.numeric(0)
      datelist<-split(v,v$Date)
      for (i in datelist){
        if((unique(i$Date)-6>=min(vac$Date))&(unique(i$Date)<=max(vac$Date))){
          name<-unique(i$vaccine)
          total<-sum(i$total_vaccinations)
          for(k in 1:length(name)){
            vac[(which(vac$Date==unique(i$Date))-6):which(vac$Date==unique(i$Date)),name[k]]<-
              as.numeric(i$total_vaccinations[which(i$vaccine==name[k])])/total
          }
        }
      }
      print(paste(unique(v$location),"has been processed"))
      return(vac)
    }
    return(NULL)
  })
  list2<-lapply(c("United Kingdom","Israel"), function(i){
    vac<-subset(V1,V1$location==i)
    vac$total_vaccinations<-cummax(vac$total_vaccinations)
    vac$Date<-as.Date(vac$Date)
    if (nrow(vac)>0){
      vac$'Johnson&Johnson'<-as.numeric(0)
      vac$'Moderna'<-as.numeric(0)
      vac$'Oxford/AstraZeneca'<-as.numeric(0)
      vac$'Pfizer/BioNTech'<-as.numeric(0)
      vac$'Sinovac'<-as.numeric(0)
      vac$'CanSino'<-as.numeric(0)
      vac$'Sputnik V'<-as.numeric(0)
      vac$'Sinopharm/Beijing'<-as.numeric(0)
      if(i=="United Kingdom"){
        vac$Moderna[which(vac$Date>="2021-04-07")]<-17/157
        vac$`Oxford/AstraZeneca`[which(vac$Date>="2021-04-07")]<-100/157
        vac$`Oxford/AstraZeneca`[which(vac$Date<"2021-04-07")]<-100/140
        vac$`Pfizer/BioNTech`[which(vac$Date>="2021-04-07")]<-40/157
        vac$`Pfizer/BioNTech`[which(vac$Date<"2021-04-07")]<-40/140
      }else{
        vac$`Pfizer/BioNTech`<-1
      }
    return(vac)
    }
  })
  data<-do.call(rbind,c(list,list2))[,-3]
  return(data)
}
pop_age_calculate<-function(agedata){
  age<-read.csv(agedata,stringsAsFactors = F)
  aa<-do.call(rbind,lapply(split(age,age$Country.or.Area),FUN=function(a){
    if(length(unique(a$Source.Year))>1){
      a<-a[!duplicated(a$Age), ]
    }
    total<-a$Value[which(a$Age=="Total")]
    a$Age<-as.numeric(a$Age)
    age0_17<-sum(a$Value[which(a$Age<18)])
    age18_60<-sum(a$Value[which(a$Age>=18&a$Age<60)])
    age60up<-sum(a$Value[which(a$Age>=60)])
    d<-as.data.frame(do.call(cbind,list(unique(a$Country.or.Area))))
    d$age0_17<-round(age0_17/total,2)
    d$age18_60<-round(age18_60/total,2)
    d$age60up<-round(age60up/total,2)
    return(d)
  }))
  return(aa)
}
RatioVariable_merge<-function(Dataset,vaccine_effect,variant_tran,list){
  #calculte the efficient vaccination rate
  effect<-read.csv(vaccine_effect,stringsAsFactors = F)
  vv<-read.csv(variant_tran,stringsAsFactors = F)
  effect[,2:ncol(effect)]<-effect[,2:ncol(effect)]/100
  effect_mean<-effect[,c("Alpha","Beta","Gamma","Delta",
                         "Eta","Kappa","base")]
  effect_sd<-as.data.frame(do.call(cbind,list(approx_sd(effect$AIQR.,effect$AIQR..1),
                                              approx_sd(effect$BIQR.,effect$BIQR..1),
                                              approx_sd(effect$GIQR.,effect$GIQR..1),
                                              approx_sd(effect$DIQR.,effect$DIQR..1),
                                              approx_sd(effect$EIQR.,effect$EIQR..1),
                                              approx_sd(effect$KIQR.,effect$KIQR..1),
                                              approx_sd(effect$IQR.,effect$IQR..1))))
  VT_mean<-vv$R0_transmiss
  VT_sd<-as.vector(approx_sd(vv$IQR.,vv$IQR..1))
  va<-as.vector(do.call(rbind,lapply(seq(1,nrow(vv)),function(k){
    rnorm(1,VT_mean[k],VT_sd[k])
  })))
  R0_base<-rnorm(length(unique(Dataset$iso_code)),3.25,0.5) #Generate R0 randomly
  Dataset<-do.call(rbind,lapply(list,FUN=function(c){
    dd<-split(Dataset,Dataset$iso_code)[[c]]
    dd$R0<-R0_base[c]*(dd[,8]*va[1]+dd[,9]*va[2]+dd[,10]*va[3]+dd[,11]*va[4]+
                         dd[,12]*va[5]+dd[,13]*va[6]+dd[,14])
    return(dd)
  }))
  V_alpha<-vector()
  V_beta<-vector()
  V_gamma<-vector()
  V_delta<-vector()
  V_eta<-vector()
  V_kappa<-vector()
  V_SARS_CoV_2<-vector()
  for(j in 1:8){
    V_alpha[j]<-rnorm(1,effect_mean[j,1],effect_sd[j,1])
    V_beta[j]<-rnorm(1,effect_mean[j,2],effect_sd[j,2]) 
    V_gamma[j]<-rnorm(1,effect_mean[j,3],effect_sd[j,3])
    V_delta[j]<-rnorm(1,effect_mean[j,4],effect_sd[j,4])
    V_eta[j]<-rnorm(1,effect_mean[j,5],effect_sd[j,5])
    V_kappa[j]<-rnorm(1,effect_mean[j,6],effect_sd[j,6])
    V_SARS_CoV_2[j]<-rnorm(1,effect_mean[j,7],effect_sd[j,7])
  }
  effect<-as.data.frame(do.call(cbind,lapply(seq(1,8),FUN=function(j){
    Dataset[,8]*V_alpha[j]+Dataset[,9]*V_beta[j]+
      Dataset[,10]*V_gamma[j]+Dataset[,11]*V_delta[j]+
      Dataset[,12]*V_eta[j]+Dataset[,13]*V_kappa[j]+
      Dataset[,14]*V_SARS_CoV_2[j]
  })))
  Dataset$Fully_vaccinated_effect<-Dataset$Fully_vaccinated_pre*
    (Dataset[,21]*effect[,1]+Dataset[,22]*effect[,2]+
       Dataset[,23]*effect[,3]+Dataset[,24]*effect[,4]+
       Dataset[,25]*effect[,5]+Dataset[,26]*effect[,6]+
       Dataset[,27]*effect[,7]+Dataset[,28]*effect[,8])
  return(Dataset)
}

select_data<-function(school_holiday){
  SH<-read.csv(school_holiday, stringsAsFactors = FALSE)
  SH<-subset(SH,SH$Year>=2019)
  SH$Date<-as.Date(SH$Date)+365
  SH<-SH[,c(1:3,13)]
  new<-SH
  new$Date<-as.Date(new$Date)+365
  SH<-do.call(rbind,list(SH,new))
  colnames(SH)<-c("iso_code","Date","CountryName","Holiday")
  return(SH)
}
COR<-function(N){
  p1<-quickcor(N, cor.test = TRUE)+geom_star(data = get_data(type = "upper" ,show.diag = FALSE))+
    geom_mark(data = get_data(type = "lower", show.diag = FALSE),size=2.5,family="T")+
    geom_abline(size=0.5,slope=-1,intercept =7)+
    scale_fill_gradient2(midpoint = 0, low = "#00498d", mid = "white", high = "#74230a",space="Lab")+
    theme(legend.position = "right",
          legend.text = element_text(color="black",size=9),
          legend.title = element_text(color="black",size=9),
          legend.key.height=unit(10,'mm'),
          legend.key.width=unit(2,'mm'),
          axis.text = element_text(color="black",size=9))
  return(p1)
  #ggsave("picture/cor.pdf",p1,width=200,height=200,units="mm",device = cairo_pdf)
}
Npidataset<-function(NPI){
  SH<-select_data(school_holiday)
  N<-read.csv(NPI,stringsAsFactors = F)
  N[is.na(N)]<-0
  N<-N[-which(N$RegionName>0),]
  N<-N[,c("CountryName","CountryCode","Date",
          "C1_School.closing","C2_Workplace.closing","C3_Cancel.public.events",
          "C4_Restrictions.on.gatherings","C5_Close.public.transport",
          "C6_Stay.at.home.requirements","C7_Restrictions.on.internal.movement",
          "C8_International.travel.controls", "H6_Facial.Coverings","StringencyIndex")]
  colnames(N)<-c("location","iso_code","Date",
                 "School closure",
                 "Workplace closure",
                 "Public events closure",
                 "Gathering restrictions",
                 "Public transport closure",
                 "Stay-at-home order",
                 "Internal movement restrictions",
                 "International travel controls",
                 "Facial covering","StringencyIndex")
  N$Date<-as.Date.character(N$Date,format="%Y%m%d")
  return(N)
}

merge_dataset_V2<-function(RT,Variants,Vaccination,NPI,POP,index,Env){
  Var<-variant_add(RT,Variants)
  N<-Npidataset(NPI)
  P<-read.csv(POP,stringsAsFactors = F)
  I<-read.csv(index,stringsAsFactors = F)
  I$aging<-as.numeric(sub("%","",I$aging))
  E<-read.csv(Env,stringsAsFactors = F)
  colnames(E)<-c("iso_code","Date","Hum","Tem")
  E$Date<-as.Date(E$Date)
  V1<-Vaccination_complementation(Vaccination)
  V<-V1[,c(1,2,3,5,6)]
  colnames(V)<-c("Date","location","iso_code","vaccinated","fully_vaccinated")
  V$Date<-as.Date(V$Date)
  Vman<-read.csv(Vaccman,stringsAsFactors = F)
  vaccine_ratio<-vaccine_ratio_cal(Vman,V1)
  Dataset<-do.call(rbind,lapply(split(Var,Var$iso_code),FUN=function(r){
    r$POP<-P$T_pop_2020[which(P$iso3c==unique(r$iso_code))]
    r$Popdensity<-I$popdensity[which(I$iso3c==unique(r$iso_code))]
    r$HealthIndex<-I$Health_index[which(I$iso3c==unique(r$iso_code))]
    r$Aging<-I$aging[which(I$iso3c==unique(r$iso_code))]
    vac_ratio<-subset(vaccine_ratio,vaccine_ratio$location==unique(r$location))[,-c(2,3)]
    npi<-subset(N,N$iso_code==unique(r$iso_code))[,-c(1:2)]
    e<-subset(E,E$iso_code==unique(r$iso_code))[,-1]
    v<-subset(V,V$iso_code==unique(r$iso_code))[,-c(2,3)]
    v[is.na(v)]<-0
    v$vaccinated<-cummax(v$vaccinated)
    v$fully_vaccinated<-cummax(v$fully_vaccinated)
    if(min(c(nrow(vac_ratio), nrow(r), nrow(npi), nrow(e),nrow(v)))>0){
      r<-merge(r,e,by="Date")
      r<-merge(r,vac_ratio,by="Date",all=T)
      r[is.na(r)]<-0
      r<-merge(r,npi,by="Date")
      r<-merge(r,v,by="Date",all=T)
      r<-r[0:which(r$Date==max(r$Date[which(r$vaccinated!=0)])),]
      r[is.na(r)]<-0
      r$Vaccinated_pre<-r$vaccinated/r$POP
      r$Fully_vaccinated_pre<-r$fully_vaccinated/r$POP
      r<-subset(r,r$iso_code!=0)
      r$Total_cases_per<-r$Total_cases/r$POP
      return(r)
    }
  }))
  return(Dataset)
}

alpha_calculate_stringency<-function(m2,DatasetV2,start=as.Date("2020-08-01"),
                                   end=as.Date("2021-07-29"),output){
  DatasetV2$StringencyIndex<-DatasetV2$StringencyIndex/100
  D<-subset(DatasetV2,DatasetV2$Date>=start&DatasetV2$Date<=end)
  X1<-as.matrix(D[,c("StringencyIndex","Fully_vaccinated_effect")])
  X2<-D[,c("HealthIndex")]
  #X2$Tem<-(X2$Tem-min(X2$Tem))/(max(X2$Tem)-min(X2$Tem))
  #X2$Hum<-(X2$Hum-min(X2$Hum))/(max(X2$Hum)-min(X2$Hum))
  #X2$Aging<-X2$Aging/max(X2$Aging)
  X2<-X2/max(X2)
  #X2$Popdensity<-X2$Popdensity/max(X2$Popdensity)
  X2<-as.matrix(X2)
  Rt<-as.vector(D[,"Rt"])
  R0<-as.vector(D[,"R0"])
  Cnum<-nrow(split(D,D$iso_code)[[1]])
  
  dataset<- list(m=length(Rt),n=length(split(D,D$iso_code)),
                 k2=ncol(X2), k1=ncol(X1),X1=X1,X2=X2,Rt=Rt,R0=R0,
                 Cnum=Cnum) 
  rstan_options(auto_write = TRUE)
  set.seed(13)
  fit<-rstan::sampling(object = m2,data=dataset,iter=1000,warmup=500,
                       chains=4,thin=1,control = 
                         list(adapt_delta = 0.95, max_treedepth =10))
  saveRDS(fit, file = output)
  #tp<-traceplot(fit,pars="alpha")
  return(fit)
}
alpha_calculate_stringency_con<-function(mcon,DatasetV2,start=as.Date("2020-10-01"),
                                         end=as.Date("2021-07-29"),output){
  DatasetV2$StringencyIndex<-DatasetV2$StringencyIndex/100
  D<-subset(DatasetV2,DatasetV2$Date>=start&DatasetV2$Date<=end)
  X1<-as.matrix(D[,c("StringencyIndex","Fully_vaccinated_effect")])
  X2<-as.matrix(D[,c("Tem")])
  X2<-(X2-min(X2))/(max(X2)-min(X2))
  Rt<-as.vector(D[,"Rt"])
  R0<-as.vector(D[,"R0"])
  dataset<-  list(m=length(Rt),k1=ncol(X1), k2=ncol(X2),
                  X2=X2,X1=X1,Rt=Rt,R0=R0)  
  rstan_options(auto_write = TRUE)
  set.seed(11)
  fit<-rstan::sampling(object = mcon,data=dataset,iter=1000,warmup=500,
                       chains=4,thin=1,control = 
                         list(adapt_delta = 0.95, max_treedepth =10))
  saveRDS(fit, file = output)
  #tp<-traceplot(fit,pars="alpha")
  return(fit)
}

contry_process_SI<-function(DatasetV2,source){
  mcon<-stan_model('stan_country.stan')
  if (dir.exists(source)==F){dir.create(source)
  }else{print("This path has been exists")}
  lapply(split(DatasetV2,DatasetV2$iso_code),FUN=function(v){
    out<-paste0(source,unique(v$location),"_withSI/")
    if (dir.exists(out)==F){dir.create(out)
    }else{print("This path has been exists")}
    start=as.Date("2020-08-01")
    fitlist<-lapply(seq(1,14),function(x){
      start<-start+30*(x-1)
      end<-start+30
      output<-paste0(out,"30days_before",end,".rds")
      fit<-alpha_calculate_stringency_con(mcon,v,start,end,output)
      return(fit)
    })
    stringencyresult(v,start,out)
  })
}

dataprocess<-function(RT,Variants,Vaccination,NPI,POP,index,Env,
                      vaccine_effect,variant_tran,google_mobility,i,outpath){
  Dataset<-merge_dataset_V2(RT,Variants,Vaccination,NPI,POP,index,Env)
  Dataset<-subset(Dataset,Dataset$continent=="Europe"|Dataset$location=="Israel")
  gc()
  DatasetV2<-RatioVariable_merge(Dataset,vaccine_effect,variant_tran,list=c(1:length(unique(Dataset$location))))
  write.csv(DatasetV2,"dataset/Europe_dataset_0923V",i,".csv",row.names = F)
  contry_process_SI(DatasetV2,source=paste0(outpath,"/Countryresult0923_withTem_V",i,"/"))
}
stringencyresult<-function(DatasetV2,start,out){
  fitlist<-list.files(out,pattern="*.rds$")
  result<-do.call(rbind,lapply(c(1:length(fitlist)),FUN = function(x){
    f<-readRDS(paste0(out,"/",fitlist[[x]]))
    start<-start+30*(x-1)
    end<-start+30
    #traceplot(f, pars = c("alpha","beta"), inc_warmup = TRUE, nrow = 4)
    k<-subset(DatasetV2,DatasetV2$Date>=start&DatasetV2$Date<=end)
    X1<-mean(k$StringencyIndex/100)
    X2<-mean(k$Fully_vaccinated_effect)
    #X3<-mean(k$Total_cases_per)
    X<-do.call(cbind,list(X1,X2))
    out1<-rstan::extract(f)
    alpha<-as.matrix(out1$alpha)
    #beta<-as.matrix(out1$beta)
    alphalist<-lapply(seq(1:ncol(alpha)),FUN = function(k){alpha[,k]*X[k]})
   # betalist<-beta[,2]*X[3]
    #alpha<-do.call(cbind,list(do.call(cbind,alphalist),betalist))
    alpha<-do.call(cbind,alphalist)
    allNPI<-apply(alpha,1,sum)
    list<-list(alpha,allNPI)
    alpha<-do.call(cbind,list)
    alpha<-(1-exp(-alpha))*100
    colnames(alpha)<-c("Stringency Index","Fully vaccinated","All NPIs")
    data<-mcmc_intervals_data(alpha,prob = .5,prob_outer= .95,point_est="median")
    data$strength<-c(X,NA)
    data$R0<-mean(k$R0)
    data$Rt<-mean(k$Rt)
    data$start<-start
    data$end<-end
    return(data)
  }))
  write.csv(result,paste0(out,"/result.csv"))
}
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
r2<-function(lm){
  pred <- predict(lm)
  n <- length(pred)
  res <- resid(lm)
  w <- weights(lm)
  if (is.null(w)) w <- rep(1, n)
  rss <- sum(w * res ^ 2)
  resp <- pred + res
  center <- weighted.mean(resp, w)
  r.df <- summary(lm)$df[2]
  int.df <- 1
  tss <- sum(w * (resp - center)^2)
  r.sq <- 1 - rss/tss
  adj.r.sq <- 1 - (1 - r.sq) * (n - int.df) / r.df
  return(r.sq)
}
expression_equationV1<-function(lm,m1){
  r.sq<-round(modelr::rsquare(lm, m1),3)
  par<-as.data.frame(summary(lm)$coefficients)
  a<-round(par$Estimate[1],2)
  eq<-substitute(italic(y)==italic(x)^a~","~~italic(r)^2~"="~r2,
                 list(a=format(a),
                      r2=format(r.sq)))
  return(eq)
}
expression_equationV2<-function(lm,m1){
  r.sq<-round(modelr::rsquare(lm, m1),3)
  mse<-
  par<-as.data.frame(summary(lm)$coefficients)
  a<-round(par$Estimate[1],2)
  b<-round(par$Estimate[2],2)
  eq<-substitute(italic(y)==a%.%italic(x)^b~","~~italic(r)^2~"="~r2,
                 list(a=format(a),
                      b=format(b),
                      r2=format(r.sq)))
  return(eq)
}
strengthcalculate<-function(lm_s,effect){
  coff<-as.data.frame(summary(lm_s)$coefficients)
  a<-coff$Estimate
  effect[which(effect<0)]<-0
  strength<-exp(log(effect)/a)
}
"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)) #利用transform函数为数据框mydata增加数据
            
            newdata <- rbind(plyr::arrange(transform(data, x = xmaxv), -y),plyr::arrange(transform(data, x = xminv), y))
            newdata_Polygon <- rbind(newdata, newdata[1,])
            newdata_Polygon$colour<-NA
            
            newdata_Path <- plyr::arrange(transform(data, x = xmaxv), -y)
            
            ggplot2:::ggname("geom_flat_violin", grobTree(
              GeomPolygon$draw_panel(newdata_Polygon, panel_scales, coord),
              GeomPath$draw_panel(newdata_Path, panel_scales, coord))
            )
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )
meta_analysis<-function(path){
  list<-list.files(path)
  conlist<-list.files(paste0(path,"/",list[1]))
  datelist<-list.files(paste0(path,"/",list[1],"/",conlist[1]),".rds")
  condata<-lapply(datelist,function(D){
    end<-substring(D,14,23)
    com<-lapply(conlist,function(c){
      result<-read.csv(paste0(path,"/",list[1],"/",c,"/result.csv"))
      X<-result$strength[which(result$end==end)][-3]
      alpha<-do.call(rbind,lapply(list,function(l){
        f<-readRDS(paste0(path,"/",l,"/",c,"/",D))
        alpha<-as.data.frame(rstan::extract(f)$alpha)
        for (k in seq(1:ncol(alpha))){
          alpha[,k]<-alpha[,k]*X[k]}
        alpha[,3]<-alpha[,1]+alpha[,2]
        return(alpha)
      }))
      d<-mcmc_intervals_data(alpha,prob = .5,prob_outer= .95,point_est="median")
      data<-data.frame(n=nrow(alpha),sd_SI=sd(alpha$V1),V_SI=sd(alpha$V2),
                       min_SI=min(alpha$V1),max_SI=max(alpha$V1),
                       min_V=min(alpha$V2),max_V=max(alpha$V2),
                       mean_SI=mean(alpha$V1),mean_V=mean(alpha$V2),
                       strength_SI=X[1],median_SI=d$m[1],q1_SI=d$ll[1],q3_SI=d$hh[1],
                       strength_V=X[2],median_V=d$m[2],q1_V=d$ll[2],q3_V=d$hh[2],
                       strength_Join=X[1],median_Join=d$m[3],q1_Join=d$ll[3],q3_Join=d$hh[3],
                       country=strsplit(c,"_")[[1]][1],
                       start=as.Date(end)-30,R0=unique(result$R0[which(result$end==end)]),
                       Rt=unique(result$Rt[which(result$end==end)]))
      return(data)
    })
    f<-do.call(rbind,com)
    if(max(f$median_SI)>0){
      ds<-subset(f,f$median_SI>0)
      sim<-metamean(n=ds$n,median=ds$median_SI,q1=ds$q1_SI,
                    q3=ds$q3_SI,studlab = ds$country)
      si<-summary(sim)$random
      s<-data.frame(start=f$start[1],par="Stringency Index",
                    strength=mean(f$strength_SI),
                    m=si$TE,se=si$seTE,lower=si$lower,upper=si$upper)}else{s<-NULL}
    if(max(f$median_V)>0){
      dv<-subset(f,f$median_V>0)
      vm<-metamean(n=dv$n,median=dv$median_V,q1=dv$q1_V,
                   q3=dv$q3_V,studlab = dv$country)
      va<-summary(vm)$random
      v<-data.frame(start=dv$start[1],par="Vaccination",
                    strength=mean(dv$strength_V),
                    m=va$TE,se=va$seTE,lower=va$lower,upper=va$upper)
    }else{v<-NULL}
    X<-do.call(rbind,list(s,v))
    return(list(f,X))
  })
  d1<-do.call(rbind,lapply(condata,function(k){as.data.frame(k[[2]])}))
  write.csv(d1,paste0(figpath,"/fig2datanew.csv"),row.names = F)
  contryresult<-do.call(rbind,lapply(condata,function(k){as.data.frame(k[[1]])}))
  write.csv(contryresult,paste0(figpath,"/countryresult_all.csv"),row.names = F)
}
