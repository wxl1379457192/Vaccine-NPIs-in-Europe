
fig1a1plot<-function(DatasetV2,figpath){
  D<-subset(DatasetV2,DatasetV2$Date>="2020-08-01")
  D<-D[,c("Date","iso_code","New_cases_smoothed")]
  D$id<-1:nrow(D)
  angle<-90-360*(D$id-0.5)/nrow(D)
  basean <- do.call(rbind,lapply(split(D,D$iso_code),FUN = function(v){
    f<-data.frame(start=min(v$id),end=max(v$id),iso_code=unique(v$iso_code))
    f$mid<-(f$start +f$ end) / 2
    return(f)
  }))
  textangle<-angle[floor(basean$mid)]
  colorpal<-c("#a18480","#5e7b7f","#a99fc6","#d2691e","#5e86c1","#a1793e",
              "#808000","#ebbe00","#8a322e","#d9e0c4","#506797","#bfac85",
              "#9084b4","#9c9c9c","#7f6358","#a0522d","#c4daee","#675f4c",
              "#c3aebf","#3c5140","#bcd1a2","#354b36","#868381","#b8977b",
              "#7c5335","#cac8b5","#fff1cd")
  g<-ggplot(data=D,aes(id, New_cases_smoothed/1000,fill=iso_code))+
    geom_col(position = position_dodge2(),width = 0.9)+
    geom_segment(data=basean,aes(x=start,y=-5,xend=end,yend=-5,color=iso_code))+
    geom_text(data = basean, aes(x = mid, y = -10, label = iso_code), 
              angle=textangle,colour = "grey40",size=2) +
    coord_polar()+
    scale_fill_manual(values=colorpal)+
    scale_color_manual(values=colorpal)+
    scale_y_continuous(limits=c(-100,60),expand=c(0,0))+
    theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          legend.position = "",
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())
  ggsave(g, filename = paste0(figpath,"/fig1a-1.pdf"),unit="mm",device = cairo_pdf)
}
fig1a2plot<-function(DatasetV2,figpath){
  world<-map_data("world")
  Dataset<-do.call(rbind,lapply(split(DatasetV2,DatasetV2$location),function(x){
    data.frame(iso=unique(x$iso_code),region=unique(x$location),Vaccination_Fully=max(x$Fully_vaccinated_pre))
  }))
  world<-full_join(world,Dataset,by="region")
  world$w<-1
  world$w[which(is.na(world$iso_code))]<-0
  world$w<-as.factor(world$w)
  p1<-ggplot()+
    geom_polygon(data = world,aes(x=long,y=lat,group=group,fill=Vaccination_Fully),color="white",size=0.1, inherit.aes = F)+  
    scale_fill_gradientn(colors =c("#bfcee5","#c5c8e8","#5b5fbf"),na.value  = "#ededed",)+
    coord_map("ortho",orientation = c(50,10,0))+
    theme(legend.position = "bottom",
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.border = element_blank())
  ggsave(p1, filename = paste0(figpath,"/fig1a-2.pdf"),width=300,heigh=300,unit="mm",device = cairo_pdf)
}
fig1bplot<-function(DatasetV2,figpath){
  cpcol1 <- c("#697fa0", "#526582", "#4682b4", "#4682b4", "#5e86c1", 
              "#97c0e2", "#5b94c2", "#80acd0", "#2f5a7f", "#0f3c62", 
              "#6587bf", "#9bb1d4", "#36517d", "#7b8ca9", "#b3bccc", 
              "#8d8fb3", "#8674a1", "#8d8fb3", "#666999", "#cacbdd", 
              "#78699f", "#b9b1cd", "#778899", "#719fda", "#496e9c",
              "#5686bf","#a099b2")
  c<-DatasetV2[,c("Date","location","StringencyIndex","Fully_vaccinated_pre","R0","Rt")]
  c$Fully_vaccinated_pre<-c$Fully_vaccinated_pre*100
  cs<-Dataset[,c("Date","location","StringencyIndex")]
  f<-ggplot()+
    geom_step(data=c,aes(x=as.Date(Date),y=Fully_vaccinated_pre,color=location),size=0.3,alpha=0.1)+
    geom_smooth(data=c,aes(x=as.Date(Date),
                           y=Fully_vaccinated_pre),color="#606bbc",fill="#b9bad8",span=0.5,size=0.9)+
    scale_x_date(limits=as.Date(c("2020-08-01",max(Dataset$Date))),
                 breaks=("2 month"),labels=date_format("%b\n%Y"),expand=c(0,0))+
    scale_y_continuous(limits=c(0,100),expand = c(0,0),breaks=seq(0,100,20))+
    labs(y ="Intensity (%)",x = NULL)+theme_light()+
    scale_color_manual(values=cpcol1)+
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          plot.title = element_text(color="black",hjust = 0,vjust=0,size = 10),
          axis.title.x= element_text(color="black",size = 10),
          axis.text.x= element_text(color="black",size = 10),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.background=element_rect(fill = "transparent",colour = NA))
  ggsave(paste0(figpath,"/fig1b-1.pdf"),f,width=100,height=80,units="mm",device = cairo_pdf)
  
  cpcol3<-c("#d9d7cb","#e2dec2","#ece5b8","#fbf1a9","#fae780",
            "#f8c43f","#c59007","#fafb94","#f8fa70","#eaeaae",
            "#d5d661","#c29b32", "#ffe6a3","#dcce8f","#d6c499",
            "#c2b54c","#d1c13d","#f0ebc1","#d8ca88","#87772c",
            "#fdeb86","#c4a603","#fbd509","#dfd4c5","#eed6b6",
            "#e0b376","#b0a586")
  f<-ggplot()+
    geom_step(data=cs,aes(x=as.Date(Date),y=StringencyIndex,color=location),size=0.3,alpha=0.1)+
    geom_smooth(data=cs,aes(x=as.Date(Date),
                            y=StringencyIndex),linetype=2,color="#dca904",fill="#feefb9",span=0.5,size=0.9)+
    scale_x_date(limits=as.Date(c("2020-08-01",max(Dataset$Date))),
                 breaks=("2 month"),labels=date_format("%b\n%Y"),expand=c(0,0))+
    scale_y_continuous(limits=c(0,100),expand = c(0,0),breaks=seq(0,100,20))+
    labs(y ="Intensity (%)",x = NULL)+theme_light()+
    scale_color_manual(values=cpcol3)+
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          plot.title = element_text(color="black",hjust = 0,vjust=0,size = 10),
          axis.title.x= element_text(color="black",size = 10),
          axis.text.x= element_text(color="black",size = 10),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.background=element_rect(fill = "transparent",colour = NA))
  ggsave("0924picture/fig1b-2.pdf",f,width=100,height=80,units="mm",device = cairo_pdf)
}
fig1cplot<-function(DatasetV2,figpath){
  Dataset<-subset(DatasetV2,DatasetV2$Date>="2020-08-01")
  var<-Dataset[,c("Date","alphaN","betaN","gammaN","deltaN",
                  "etaN","kappaN","baseN")]
  var_mean<-group_by(var,Date)%>%summarize_each(funs(mean))
  v<-melt(var_mean,id=c("Date"))
  v1<-subset(v,v$Date=="2021-01-31")
  v1$Date<-"2021-01-30"
  v<-do.call(rbind,list(v,v1))
  v$value[which(v$value==0)]<-NA
  var<-ggplot(v,aes(x=as.Date(Date),y=variable))+
    geom_tile(aes(fill=value),
              color="white",lwd=0.1,
              linetype=1)+
    scale_fill_gradientn(colours= brewer.pal(7,"YlGnBu"),na.value="grey99",
                         values=c(0,0.01,0.3,0.6,1.0))+
    scale_x_date(limits=as.Date(c("2020-08-01",max(v$Date))),
                 breaks=("1 month"),labels=date_format("%b\n%Y"),expand=c(0,0))+
    labs(y =NULL,x = NULL)+theme_light()+
    theme(legend.position = "right",
          legend.title=element_blank(),
          legend.key.height=unit(0.4,'cm'),
          legend.key.width =unit(0.1,'cm'),
          #panel.grid=element_line(color="transparent",size=0.01),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(color="black",size = 10),
          axis.title.x= element_text(color="black",size = 10),
          axis.text.x= element_text(color="black",size = 10),
          plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
          panel.background=element_rect(fill = "transparent",colour = NA))
  ggsave(paste0(figpath,"/variant_distribution.pdf"),var,width=180,height=40,units="mm",device = cairo_pdf)
  dev.off()
}
fig1dplot<-function(DatasetV2,figpath){
  vv<-DatasetV2[,c(1,21:24,27:28)]
  vv1<-group_by(vv,Date)%>%summarize_each(funs(mean))
  vc<-melt(vv1,id=c("Date"))
  vc1<-subset(vc,vc$Date=="2021-01-31")
  vc1$Date<-"2021-01-30"
  vc<-do.call(rbind,list(vc,vc1))
  vc$value[which(vc$value==0)]<-NA
  vacman<-ggplot(vc,aes(x=as.Date(Date),y=variable))+
    geom_tile(aes(fill=value),
              color="white",lwd=0.1,
              linetype=1)+
    scale_fill_gradientn(colours= brewer.pal(7,"YlGnBu"),
                         na.value="grey99",values=c(0,0.05,0.2,0.5,0.75,1.0))+
    scale_x_date(limits=as.Date(c("2020-08-01",max(v$Date))),
                 breaks=("1 month"),labels=date_format("%b\n%Y"),expand=c(0,0))+
    labs(y =NULL,x = NULL)+theme_light()+
    theme(legend.position = "right",
          legend.title=element_blank(),
          legend.key.height=unit(0.3,'cm'),
          legend.key.width =unit(0.1,'cm'),
          plot.title = element_text(color="black",size = 10),
          axis.title.x= element_text(color="black",size = 10),
          axis.text.x= element_text(color="black",size = 10),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.background=element_rect(fill = "transparent",colour = NA))
  ggsave(paste0(figpath,"/vaccine_distribution.pdf"),vacman,width=180,height=35,units="mm",device = cairo_pdf)
  
}



fig2a1plot<-function(figpath){
  d1<-read.csv(paste0(figpath,"/fig2datanew.csv"),row.names = F)
  g<-ggplot(data=d1)+
    geom_xspline(aes(x=as.Date(start),y=(1-exp(-m))*100,color=par))+
    theme_clean()+
    geom_errorbar(mapping = aes_(ymin =~(1-exp(-lower))*100, ymax=~(1-exp(-upper))*100, x=~as.Date(start),color= ~par),
                  show.legend = T,alpha=0.9,width=10,size=0.2)+
    geom_point(aes(x=as.Date(start),y=(1-exp(-m))*100,color=par),size=0.5,shape=12)+
    scale_x_date(limits=as.Date(c("2020-07-15","2021-09-16")),
                 breaks=("2 month"),labels=date_format("%b\n%Y"),expand=c(0,0))+
    scale_color_manual(values=c("#7f97bd","#723f39"))+
    scale_fill_manual(values=c("#7f97bd","#723f39"))+
    scale_size_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
    scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100),expand=c(0,0))+
    labs(y =paste0("\u0394","Rt","(%)"),x = NULL)+
    guides(color=guide_legend(nrow=1,byrow=TRUE))+
    theme(legend.position = "bottom",
          panel.grid.major = element_line(colour = "transparent"),
          panel.grid.minor = element_line(colour = "transparent"),
          legend.title=element_blank(),
          plot.title = element_text(color="black",hjust = 0.5,vjust=0,size = 10),
          axis.title.x= element_blank(),
          axis.text.x = element_text(color="black",hjust=0.5,size = 10),
          axis.title.y= element_text(color="black",size = 10),
          axis.text.y= element_text(color="black",hjust=1,size = 10),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.background=element_rect(fill = "transparent",colour = NA))
  ggsave(paste0(figpath,"/fig2a.pdf"),g,units="mm",width=90,height=90,device = cairo_pdf)
  dev.off()
}
fig2a2plot<-function(DatasetV2,figpath){
  Dataset<-subset(DatasetV2,DatasetV2$Date>="2020-08-01")
  c<-Dataset[,c("Date","location","Rt")]
  R<-Dataset[,c("Date","R0","Rt")]
  R1<-R%>%dplyr::group_by(Date)%>%summarize_each(funs(mean))
  cbPalette <- c("#f5ead6", "#e6c68e", "#deb368", "#d1952e", "#8e661f", 
                 "#d3c6b1", "#bba886", "#a78e62", "#89734d", "#52452e", 
                 "#d8d1ca", "#c7bdb3", "#aa9b8d", "#7e6d5d", "#726040", 
                 "#e5be9e", "#daa172", "#ce8346", "#ad672e", "#895224", 
                 "#feeab4", "#fdde87", "#fcd055", "#ebae05", "#7d5d02","#a27c3c")
  f<-ggplot(data=R1)+
    geom_smooth(aes(x=as.Date(Date),y=Rt),color="#60868d",span=0.4,size=0.9,alpha=0.4)+
    geom_smooth(aes(x=as.Date(Date),y=R0),color="#9cbcbc",span=0.4,size=0.9,alpha=0.4,linetype=4)
  gg1 <- ggplot_build(f)
  # extract data for the loess lines from the 'data' slot
  df2 <- data.frame(x = gg1$data[[1]]$x,
                    ymin = gg1$data[[1]]$y,
                    ymax = gg1$data[[2]]$y) 
  f2<-f+
    geom_ribbon(data=df2,aes(x=as.Date(x), ymin = ymin, ymax = ymax),fill="#c6e7e6",
                alpha = 0.3)+
    scale_x_date(limits=as.Date(c("2020-08-01",max(Dataset$Date))),
                 breaks=("2 month"),labels=date_format("%b\n%Y"))+
    scale_y_continuous(limits=c(0,5),expand = c(0,0),breaks=seq(0,5,1))+
    geom_vline(aes(xintercept=as.Date("2020-12-22"),color="#d2dde4",alpha=0.3))+
    geom_vline(aes(xintercept=as.Date("2021-05-04"),color="#d2dde4",alpha=0.3))+
    labs(y ="Repoduction number",x = NULL)+theme_clean()+
    scale_color_manual(values=cbPalette)+
    theme(legend.position = "",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          plot.title = element_text(color="black",hjust = 0,vjust=0,size = 10),
          axis.title.x= element_text(color="black",size = 10),
          axis.text.x= element_text(color="black",size = 10),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.background=element_rect(fill = "transparent",colour = NA))
  ggsave(paste0(figpath,"/fig2a-2.pdf"),f2,width=90,height=30,units="mm",device = cairo_pdf)
}
fig3plot<-function(figpath){
  rd<-read.csv(paste0(figpath,"/countryresult_all.csv"),stringsAsFactors = F)
  l<-rd[,c("median_SI","strength_SI","median_V","strength_V","start","country")]
  l1<-rd[,c("median_SI","strength_SI","strength_V")]
  colnames(l1)<-c("m","s","v")
  l2<-rd[,c("median_V","strength_SI","strength_V")]
  colnames(l2)<-c("m","s","v")
  l1$par<-"Stringency Index"
  l2$par<-"Fully vaccination"
  l<-do.call(rbind,list(l1,l2))
  l$v1<-cut(l$v*100,breaks=seq(0,60,15))
  l$s1<-cut(l$s*100,breaks=seq(20,80,15))
  l<-l[complete.cases(l),]
  l$m<-(1-exp(-l$m))
  l$label<-"After vaccination"
  l$label[which(l$v==0)]<-"Before vaccination"
  lfig<-subset(l,l$v>0)
  lfig3<-subset(l,l$v==0)
  lfig3<-subset(lfig3,lfig3$par=="Stringency Index")
  lfig1<-subset(lfig,lfig$par=="Stringency Index")
  lfig1<-subset(lfig1,lfig1$m>0.01)
  lfig2<-subset(lfig,lfig$par=="Fully vaccination")
  my_comparisons <- list(c("(0,15]", "(15,30]"), c("(0,15]", "(30,45]"), c("(0,15]", "(45,60]"),
                         c("(15,30]", "(30,45]"), c("(15,30]", "(45,60]"), c("(30,45]", "(45,60]")) 
  #my_comparisons <- list(c("(0,20]", "(20,40]"), c("(0,20]", "(40,60]"),  c("(20,40]", "(40,60]"))
  #my_comparisons <- list(c("(0,10]", "(10,20]"), c("(0,10]", "(20,30]"), c("(0,10]", "(30,40]"),c("(0,10]", "(40,50]"), c("(0,10]", "(50,60]"),
  #                       c("(10,20]", "(20,30]"), c("(10,20]", "(30,40]"),c("(10,20]", "(40,50]"), c("(10,20]", "(50,60]"),
  #                       c("(20,30]", "(30,40]"), c("(20,30]", "(40,50]"), c("(20,30]", "(50,60]"),
  #                       c("(30,40]", "(40,50]"),c("(30,40]", "(50,60]"),c("(40,50]", "(50,60]")) 
  fig3a<-ggplot()+
    geom_flat_violin(data=lfig3,aes(x=as.factor(v1),y=m*100),fill ="#7f97bd", 
                     alpha=0.4, color="black",width=1)+
    geom_flat_violin(data=lfig1,aes(x=as.factor(v1),y=m*100),fill ="#7f97bd", 
                     alpha=0.4, color="black",width=1)+
    geom_boxplot(data=lfig2,aes(x=as.factor(v1),y=m*100),width=0.2,position=position_nudge(x=0.5),color="#723f39",fill="#d4afaa",alpha=0.4)+
    geom_jitter(data=lfig1,aes(x=as.factor(v1),y=m*100,color=s),position=position_nudge(x=-0.05),size=0.5,alpha=0.8)+
    geom_jitter(data=lfig3,aes(x=as.factor(v1),y=m*100,color=s),position=position_nudge(x=-0.05),size=0.5,alpha=0.8)+
    scale_color_gradientn(colours=brewer.pal(7,"PuBuGn"),breaks=c(min(lfig1$s),max(lfig1$s)))+
    scale_fill_gradientn(colours=  c("grey95","#7b2625"),breaks=c(0,max(lfig2$v)))+
    labs(y = paste0("Reduction in R0,t","(%)"),x = "Vaccination rate (%)")+
    scale_y_continuous(limits=c(0,80),expand=c(0,0))+
    theme_clean()+
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          plot.title = element_text(color="black",hjust = 0.5,vjust=0,size = 10),
          axis.title.x= element_text(color="black",size = 10),
          axis.text.x = element_text(color="black",size = 10),
          axis.title.y= element_text(color="black",size = 10),
          axis.text.y = element_text(color="black",size = 10),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.background=element_rect(fill = "transparent",colour = NA))
  ggsave(paste0(figpath,"/fig3a.pdf"),fig3a,units="mm",width=120,height=100,device = cairo_pdf)
  
  mydata2<-do.call(rbind,lapply(split(lfig1,lfig1$s1),function(data){
    k<-do.call(rbind,lapply(my_comparisons,function(m){
      d<-do.call(rbind,split(data,data$v1)[m])
      if(length(unique(d$v1))>1){
        c<-compare_means( m ~ v1, data = d
                          #,method="anova" 
                          )
        #c<-anova_test(m ~ v1,data=d)
        c$par<-paste(as.vector(m)[1],as.vector(m)[2])
      }else{c<-NULL}
      return(c)
    }))
    return(k)
  }))
  write.csv(mydata2,paste0(figpath,"/fig3b.csv"),row.names = F)
  fig3b<-ggplot()+
    geom_boxplot(data=lfig1,aes(x=s1,y=m*100,fill=v1,color=v1), 
                 alpha=0.5,width=1)+
    scale_fill_manual(values = c("#cec2d6","#fde590",
                                 "#f7a95f","#14a1a7",
                                 "#9bbeda","#6e6e4b",
                                 "#e3c6c2","#c4c3be"
    ))+
    scale_color_manual(values = c("#cec2d6","#fde590",
                                  "#f7a95f","#14a1a7",
                                  "#9bbeda","#6e6e4b",
                                  "#e3c6c2","#c4c3be"))+
    labs(y = paste0("Reduction in R0,t","(%)"),
         x = "Stringency index intensity")+
    theme_clean()+
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          plot.title = element_text(color="black",hjust = 0.5,vjust=0,size = 10),
          axis.title.x= element_text(color="black",size = 10),
          axis.text.x = element_text(color="black",size = 10),
          axis.title.y= element_text(color="black",size = 10),
          axis.text.y = element_text(color="black",size = 10),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.background=element_rect(fill = "transparent",
                                        colour = NA))
  ggsave(paste0(figpath,"/fig3b.pdf"),fig3b,units="mm",width=100,height=120,device = cairo_pdf)
  
  
  ################plot fig3a-2##############
  D<-DatasetV2[,c("Fully_vaccinated_effect","Date","location")]
  D$Date<-as.Date(D$Date)
  D$m<-cut(D$Date,breaks=seq(as.Date("2020-08-01"),max(D$Date)+30,30))
  D<-D[complete.cases(D),]
  D<-D%>%dplyr::group_by(m,location)%>%dplyr::summarise(Fully_vaccinated_effect=mean(Fully_vaccinated_effect))
  D$v1<-cut(D$Fully_vaccinated_effect*100,breaks=seq(-10,90,10))
  
  D$m<-as.Date(D$m)
  conv<-D[,c("m","v1","location")]%>%dplyr::group_by(v1,location)%>%dplyr::summarise(min=min(m),max=max(m))
  conv<-conv[complete.cases(conv),]
  Variants[is.na(Variants)]<-0
  Variants$date<-as.Date(Variants$date)
  Variants<-do.call(rbind,split(Variants,Variants$location)[unique(D$location)])
  
  dv<-lapply(seq(1,10),function(n){
    v<-split(conv,conv$v1)[[n]]
    if(n>1){
      v1<-split(conv,conv$v1)[[n-1]]
    }
    if (nrow(v)>0){
      vlist<-lapply(1:length(split(v,v$location)),function(m){
        l<-split(v,v$location)[[m]]
        if(n>1){min<-v1$max[which(v1$location==l$location)]}else(min<-l$min)
        k<-subset(Variants,Variants$date<=l$max&Variants$date>l$min)
        k<-k[,-c(1,2)]%>%dplyr::group_by(variant)%>%dplyr::summarise_each(funs(sum))
        return(k)
      })
      i<-do.call(rbind,vlist)
      if (nrow(i>0)){
        i<-i[,-3]%>%dplyr::group_by(variant)%>%dplyr::summarise_each(funs(sum))
        alphaN<- i$num_sequences[which(i$variant=="Alpha")]/i$num_sequences_total[1]
        betaN<-i$num_sequences[which(i$variant=="Beta")]/i$num_sequences_total[1]
        gammaN<-i$num_sequences[which(i$variant=="Gamma")]/i$num_sequences_total[1]
        deltaN<-i$num_sequences[which(i$variant=="Delta")]/i$num_sequences_total[1]
        etaN<-i$num_sequences[which(i$variant=="Eta")]/i$num_sequences_total[1]
        kappaN<-i$num_sequences[which(i$variant=="Kappa")]/i$num_sequences_total[1]
        baseN<-1-alphaN-betaN-gammaN-deltaN-etaN-kappaN
        r<-data.frame(alpha=alphaN,beta=betaN,
                      gamma=gammaN,delta=deltaN,eta=etaN,kappa=kappaN,
                      base=baseN,v=as.vector(unique(v$v1)))
      }else{r<-NULL}
    }else{r<-NULL}
    return(r)
  })
  varD<-do.call(rbind,dv)
  varD$v<-unique(conv$v1)[-c(8,9)]
  varD<-melt(varD,id="v")
  fig2b4<-ggplot(data=varD,aes(v,value,fill=variable))+geom_bar(stat="identity",position="stack")+
    scale_fill_manual(values=  c("#d2dde4","#f5ba93","#80bfb0","#ead4e2","#75cbe8","#f6dd6e","#fff7ee"))+
    labs(y = paste0("Proportion of variant","(*100%)"),x = "Vaccination rate (%)")+
    scale_y_continuous(limits=c(0,1),expand=c(0,0))+
    theme_clean()+
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          plot.title = element_text(color="black",hjust = 0.5,vjust=0,size = 10),
          axis.title.x= element_text(color="black",size = 10),
          axis.text.x = element_text(color="black",size = 10),
          axis.title.y= element_text(color="black",size = 10),
          axis.text.y = element_text(color="black",size = 10),
          plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"),
          panel.background=element_rect(fill = "transparent",colour = NA))
  ggsave(paste0(figpath,"/fig3a-2.pdf"),fig2b4,units="mm",width=120,height=50,device = cairo_pdf)

}
fig4plot<-function(risk,figpath){
  risk<-do.call(rbind,split(risk,risk$CountryName)[unique(DatasetV2$location)])
  risk$Date<-as.Date(risk$Date)
  risk<-subset(risk,risk$Date>="2020-08-01")
  rd<-read.csv(paste0(figpath,"/countryresult_all.csv"),stringsAsFactors = F)
  l<-rd[,c("median_SI","strength_SI","median_V","strength_V","start","country")]
  dda<-max(DatasetV2$Date)
  Rtdata<-DatasetV2[,c("iso_code","Date","R0","Rt","Fully_vaccinated_effect","StringencyIndex")]
  Rtdata<-subset(Rtdata,Rtdata$Date==dda)
  Rtdata<-Rtdata%>%distinct(iso_code,Date,.keep_all = T)
  Rt<-Rtdata
  Rt<-Rtdata[,-2]%>%dplyr::group_by(iso_code)%>%dplyr::summarise_each(funs(mean))
  Rt$Reduction<-(Rt$Rt-1)/Rt$R0*100
  Rt$Open<-abs(Rt$Reduction)^(1/0.93)
  Rt$Open[which(Rt$Reduction<0)]<--Rt$Open[which(Rt$Reduction<0)]
  write.csv(Rt,paste0(figpath,"/fig4a_result.csv"),row.names = F)
  colnames(Rt)<-c("adm0_a3",names(Rt)[-1])
  world <- ne_countries(scale = "medium",type = 'map_units', returnclass = "sf")
  class(world)
  world<-full_join(world,Rt,by="adm0_a3")
  world<-world[-which(world$adm0_a3=="ATA"),]
  p1<-ggplot() +geom_sf(data = world,aes(fill=Open), inherit.aes = F)+  
    scale_fill_gradientn(colors =c("#51a1b5","#fddc8c","#a62b1f"),na.value  = "#ededed")+
    coord_sf(expand = FALSE,crs=4326) +
    scale_x_continuous(limits=c(-25,45))+
    scale_y_continuous(limits=c(25,75))+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill="transparent",color=NA),
          legend.key=element_rect(fill=NA),
          plot.title = element_text(color="black",hjust = 0,vjust=0, size=10),
          axis.ticks= element_line(color="black"),
          axis.text = element_text(color="black",vjust=0,hjust=0, size=9),
          panel.background=element_rect(fill = "transparent",colour = NA),
          axis.line = element_line(size=0.05),
          axis.title = element_text(color="black", size=9))
  ggsave(p1, filename = paste0(figpath,"/fig4a.pdf"),width=120,heigh=90,unit="mm",device = cairo_pdf)
  
  dda<-as.Date("2021-03-04")
  Rtdata<-DatasetV2[,c("location","Date","R0","Rt","Fully_vaccinated_effect","StringencyIndex")]
  Rtdata<-subset(Rtdata,Rtdata$Date==dda)
  Rtdata<-Rtdata%>%distinct(location,Date,.keep_all = T)
  Rt<-Rtdata
  Rt<-Rtdata[,-2]%>%dplyr::group_by(location)%>%dplyr::summarise_each(funs(mean))
  Rt$Reduction<-(Rt$Rt-1)/Rt$R0*100
  Rt$Open<-abs(Rt$Reduction)^(1/0.93)
  Rt$Open[which(Rt$Reduction<0)]<--Rt$Open[which(Rt$Reduction<0)]
  r<-subset(risk,risk$Date==dda)[,-c(1,2)]
  r1<-r[,-2]%>%dplyr::group_by(CountryName)%>%dplyr::summarise_each(funs(mean))
  colnames(r1)<-c("location",names(r1)[-1])
  r<-merge(Rt,r1,by="location")
  g2<-ggplot(r,aes(x=Open,y=openness_risk))+geom_point(aes(color=StringencyIndex,size=Fully_vaccinated_effect))+
    geom_vline(aes(xintercept=0),linetype=2)+geom_hline(aes(yintercept=0.5),linetype=2)+
    geom_text(aes(label=location), nudge_y = 0.04,size=3)+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill="transparent",color=NA),
          legend.key=element_rect(fill=NA),
          plot.title = element_text(color="black",hjust = 0,vjust=0, size=10),
          axis.ticks= element_line(color="black"),
          axis.text = element_text(color="black",vjust=0,hjust=0, size=9),
          panel.background=element_rect(fill = "transparent",colour = NA),
          axis.line = element_line(size=0.05),
          axis.title = element_text(color="black", size=9))
  ggsave(g2, filename = paste0(figpath,"/fig4b.pdf"),width=110,heigh=90,unit="mm",device = cairo_pdf)
  return(r)
}
McmcParplot<-function(path){
  allfit<-readRDS(paste0(path,"/Countryresult0923_withTem_V1/Austria_withSI/30days_before2020-08-31.rds"))
  out_rhat<-as.data.frame(rhat(allfit))#Rhat
  colnames(out_rhat)<-"x"
  Rneff<-as.data.frame(neff_ratio(allfit))####ESS
  colnames(Rneff)<-"x"
  efdata<-list()
  list<-list.files(path)
  conlist<-list.files(paste0(path,list[1]))
  datelist<-list.files(paste0(path,list[1],"/",conlist[1]),".rds")
  ValD<-lapply(conlist,function(c){
    com<-lapply(datelist,function(D){
      dateresult<-do.call(rbind,lapply(list,function(l){
        f<-readRDS(paste0(path,l,"/",c,"/",D))
        rhat<-as.data.frame(rhat(f))
        Rneff<-as.data.frame(neff_ratio(f))
        dd<-do.call(cbind,list(rhat,Rneff))
        colnames(dd)<-c("rhat","R_ESS")
        return(dd)
      }))
      return(dateresult)
    })
    c1<-do.call(rbind,com)
    return(c1)
  })
  c<-do.call(rbind,condata)
  valplot<-list()
  valplot[[1]]<-ggplot(ValD,aes(rhat))+
    geom_histogram(bins=100,colour="#a23434",fill="#a23434")+
    labs(x=NULL,y=NULL,title="Rhat")+
    scale_x_continuous(expand=(c(0,0)),limits=c(0.995,1.060))+
    scale_y_continuous(expand=(c(0,0)))+geom_vline(xintercept = 1,linetype ="dotdash",size=1)+
    theme(axis.line.y= element_blank(),
          axis.text.y = element_blank(),
          axis.text.x =element_text(size = 8),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.title=element_text(size = 10,hjust=0),
          panel.grid=element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.ticks.x = element_line(colour = "black"),
          panel.background=element_rect(fill = "transparent",colour = NA),
          plot.margin=unit(c(0.1,0.2,0.1,0.2),"cm"),
          axis.text=element_text(size = 10))
  valplot[[2]]<-ggplot(ValD,aes(R_ESS))+
    geom_histogram(bins=100,colour="#a23434",fill="#a23434")+
    labs(x=NULL,y=NULL,title="R_ESS")+scale_x_continuous(expand=(c(0,0)),limits=c(0,1.1))+
    scale_y_continuous(expand=(c(0,0)))+geom_vline(xintercept = 1,linetype ="dotdash",size=1)+
    theme(axis.line.y= element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          axis.text.x =element_text(size = 8),
          plot.title=element_text(size = 10,family="Times New Roman",hjust=0),
          panel.grid=element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.ticks.x = element_line(colour = "black"),
          panel.background=element_rect(fill = "transparent",colour = NA),
          plot.margin=unit(c(0.1,0.2,0.1,0.2),"cm"),
          axis.text=element_text(size = 10,family="Times New Roman"))
  ggsave(paste0(figpath,"/MCMC.pdf"),grid.arrange(arrangeGrob(grobs = valplot, ncol = 1)),
         width=180,height=60,units="mm",device = cairo_pdf)
}
country_result_plot<-function(path,figpath){
  if (dir.exists(figpath)==F){dir.create(figpath)
  }else{print("This path has been exists")}
  efdata<-list()
  list<-list.files(path)
  conlist<-list.files(paste0(path,"/",list[1]))
  datelist<-list.files(paste0(path,"/",list[1],"/",conlist[1]),".rds")
  condata<-lapply(conlist,function(c){
    com<-lapply(datelist,function(D){
      dateresult<-do.call(rbind,lapply(list,function(l){
        result<-read.csv(paste0(path,"/",l,"/",c,"/result.csv"))
        end<-substring(D,14,23)
        f<-readRDS(paste0(path,"/",l,"/",c,"/",D))
        alpha<-as.data.frame(rstan::extract(f)$alpha)
        X<-result$strength[which(result$end==end)][-3]
        for (k in seq(1:ncol(alpha))){
          alpha[,k]<-alpha[,k]*X[k]}
        return(alpha)
      }))
      colnames(dateresult)<-c("Stringency index","Vaccination")
      d<-mcmc_intervals_data(dateresult,prob = .5,prob_outer= .95,point_est="median")
      d$country<-strsplit(c,"_")[[1]][1]
      d$start<-substring(D,14,23)
      s1<-read.csv(paste0(path,"/",list[[1]],"/",c,"/result.csv"))
      d$strength<-s1$strength[which(s1$start==as.Date(substring(D,14,23))-30)][-3]
      return(d)
    })
    c1<-do.call(rbind,com)
    return(c1)
  })
  c<-do.call(rbind,condata)
  write.csv(c,paste0(figpath,"/fig4data.csv"),row.names = F)
  condata<-split(c,c$country)
  plot<-lapply(condata,function(r){
    g<-ggplot(data=r)+
      geom_xspline(aes(x=as.Date(start),y=(1-exp(-m))*100,color=parameter))+
      geom_col(aes(x=as.Date(start),y=strength*100,fill=parameter ),alpha=0.1)+
      theme_clean()+
      geom_errorbar(mapping = aes_(ymin =~(1-exp(-l))*100, ymax=~(1-exp(-h))*100, x=~as.Date(start),color= ~parameter),
                    show.legend = T,alpha=0.9,width=10,size=0.2,position =  position_dodge(width=0.2))+
      geom_point(aes(x=as.Date(start),y=(1-exp(-m))*100,color=parameter ),size=0.5,shape=12)+
      #geom_smooth(aes(x=as.Date(start),y=real),span=0.4,color="#8d8e96",linetype = "dashed")+
      scale_x_date(limits=as.Date(c("2020-07-15","2021-09-20")),
                   breaks=("4 month"),labels=date_format("%b\n%Y"),expand=c(0,0))+
      scale_color_manual(values=c("#7f97bd","#723f39"))+
      scale_fill_manual(values=c("#7f97bd","#723f39"))+
      scale_size_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
      scale_y_continuous(limits=c(0,145),breaks=c(0,25,50,75,100),expand=c(0,0))+
      labs(y =paste0("Reduction in Rt","(%)"),x = NULL,title=unique(r$country))+
      guides(color=guide_legend(nrow=1,byrow=TRUE))+
      theme(legend.position = "",
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.title = element_blank(),
            plot.title = element_text(color="black",hjust = 0.5,vjust=0,size = 10,family="Times New Roman"),
            axis.title.x = element_text(color="black",size = 10),
            axis.text.x = element_text(color="black",size = 10),
            axis.ticks.y = element_line(color="black"),
            axis.title.y = element_text(color="black",angle=90,size = 10),
            axis.text.y = element_text(color="black",hjust=1,size = 10),
            plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
            panel.background=element_rect(fill = "transparent",colour = NA))
    return(g)
  })
  legend<-plot[[1]]+theme(legend.position = "bottom",
                          legend.title=element_blank(),
                          legend.text = element_text(size = 8),
                          legend.key.height=unit(0.2,'cm'),
                          legend.key.width=unit(1,'cm'),
                          legend.background = element_rect(fill="transparent"),
                          legend.key=element_rect(fill="transparent"))
  mylegend<-g_legend(legend)
  plot<-grid.arrange(mylegend,arrangeGrob(grobs =plot,ncol =2),heights = c(0.5,10))
  ggsave(paste0(figpath,"/country_result.pdf"),plot,units="mm",width=240,height=320,device = cairo_pdf)
}



