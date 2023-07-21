# global library
library(tidyverse)
library(xgboost)
library(ggplot2)
library(gridExtra)
library(h2o)

set.seed(2023)

load('crd_data/85+_weekly_mortality_CRD.rd')
load('crd_data/0-20_weekly_mortality_CRD.rd')
load('crd_data/20-40_weekly_mortality_CRD.rd')
load('crd_data/40-65_weekly_mortality_CRD.rd')
load('crd_data/65-85_weekly_mortality_CRD.rd')
load('crd_data/all_age_weekly_mortality_CRD.rd')

variable_list1 <- c("AH1", "AH3","B", "o3","no2","pm10")
variable_list2 <- c("o3","no2","pm10")
variable_list3 <- c("AH1","AH3","B")

variable_list4 <- c("AH1", "AH3","ad","rsv", "B", "o3","no2","pm10")
variable_list5 <- c("o3","no2","pm10")
variable_list6 <- c("AH1","AH3","B","ad","rsv")


df_list<-c('20','40','65','85','max','all')
df_name<-c('0-19','20-39','40-64','65-84','85+','all')

param1 <- list(
  "objective"  = "count:poisson"
  , "eval_metric" = "poisson-nloglik"
  , "eta" = .05
  , "subsample" = 0.5
  , "colsample_bytree" = 1
  , "min_child_weight" = 1
  , "max_depth" = 4
  , 'monotone_constraints'=c(0,0,0,1,1,1,1,1,1,1)
  # "week""mean_temp" "humidity" "year" "A_H1" "A-H3" "B"  "o3" "no2" "pm10" "ad","rsv  
)

plot_df1 <- data.frame(week_number = numeric(),
                       actual_death = numeric(),
                       predicted_gam = numeric(),
                       predicted_xg = numeric(),
                       age_group = character())

vaxgb1 <- data.frame(week=numeric(),
                     year=numeric(),
                     pop=numeric(),
                     variable=character(),
                     excess_death=numeric(),
                     age_group=character())

vagam1<-data.frame(week=numeric(),
                   year=numeric(),
                   pop=numeric(),
                   variable=character(),
                   excess_death=numeric(),
                   age_group=character())

paxgb1 <- list()

pagam1 <- list()


plot_imp1 <- list()


h2o.init(nthreads = -1) 

for (i in df_list){
  df <- get(paste('df',i,sep='_'))
  colnames(df)[2] <-'weekly_death'
  
  df <- df %>% 
    select(idx,weekly_death,pop,week,mean_temp,humidity,year,AH1,AH3,B,o3,no2,pm10) 
  
  df$week <- as.numeric(df$week)
  df$pop <- lowess(df$idx, df$pop)$y
  
  df <- subset(df,year>=2014)
  df_train <- subset(df,year<2020)
  
  # XGBoost model1 - AH1,AH3,B,o3,no2,pm10
  xgtrain1 <- xgb.DMatrix(as.matrix(df_train[,-c(1,2,3)]), # delete weekly_death(y),idx, pop
                          label = df_train$weekly_death)
  setinfo(xgtrain1, "base_margin", log(df_train$pop))
  xg1 <- xgb.DMatrix(as.matrix(df[,-c(1,2,3)]),
                     label = df$weekly_death)
  setinfo(xg1, "base_margin", log(df$pop))
  
  xgb1 <- xgb.cv(data = xgtrain1,
                 nrounds = 500, params = param1,nfold = 5)
  nround1 <- xgb1$evaluation_log$iter[which.min(xgb1$evaluation_log$test_poisson_nloglik_mean)]
  
  final1 <- xgboost(data = xgtrain1, params = param1, nrounds = nround1,verbose=0)
  pred_xg1 <- predict(final1,xg1)
  
  plot_imp1[[i]] <- xgb.importance(feature_names = colnames(xgtrain1),model=final1)
  
  
  
  # GAM model1 - AH1,AH3,B,ad,rsv,o3,no2,pm10
  df_h2o1 <- as.h2o(df)
  df_h2o1$death_rate <- df_h2o1$weekly_death/df_h2o1$pop # y is death rate rather than death count
  df_train1 <- df_h2o1[df_h2o1$year<2020,]
  fit1 <- h2o.gam(x=c('year','week',"AH1", "AH3", "B","no2","o3","pm10","mean_temp","humidity"),
                  y='death_rate',
                  training_frame = df_train1,
                  nfolds = 5,seed=2023,alpha = 0,family = 'gaussian',
                  gam_columns = c('week','mean_temp','humidity'),bs=c(0,0,0))
  pred1 <- h2o.predict(object = fit1, newdata = df_h2o1)
  pred_gam1 <- as.data.frame(pred1)
  
  
  plot_temp1 <- cbind(df[,c('idx','weekly_death')],pred_gam1*df$pop,pred_xg1) # get prodicted death count
  colnames(plot_temp1) <- c('week_number','actual_death','predicted_gam','predicted_xg')
  plot_temp1$age_group <- df_name[which(df_list==i)]
  plot_df1 <- rbind(plot_df1,plot_temp1)
  
  
  #variable attribute xgboost 1
  for (j in variable_list1){
    df_new <- df
    df_new[,j] <- 0
    xg <- xgb.DMatrix(as.matrix(df_new[,-c(1:3)]) 
                      # label = df$weekly_death
    )
    setinfo(xg, "base_margin", log(df_new$pop))
    pred_var <- predict(final1,xg)
    
    df_pred <- cbind(df[,c('week','year','pop')],(pred_xg1-pred_var))
    colnames(df_pred) <- c('week','year','pop','death_attr')
    df_pred$age_group <- df_name[which(df_list==i)]
    df_pred$variable <- j
    vaxgb1 <- rbind(vaxgb1,df_pred)
  }
  
  #variable attribute gam 1
  for (j in variable_list1){
    df_new <- df
    df_new[,j] <- 0
    
    df_h2o <- as.h2o(df_new)
    df_h2o$death_rate <- df_h2o$weekly_death/df_h2o$pop # y is death rate rather than death count
    
    pred_var <- h2o.predict(object = fit1, newdata = df_h2o)
    pred_var <- as.data.frame(pred_var)
    df_pred <- cbind(df[,c('week','year','pop')],((pred_gam1-pred_var)*df_new$pop))
    colnames(df_pred) <- c('week','year','pop','death_attr')
    df_pred$age_group <- df_name[which(df_list==i)]
    df_pred$variable <- j
    vagam1 <- rbind(vagam1,df_pred)
  }
  
}

write.csv(plot_df1,"model1.csv",row.names = F)

write.csv(vaxgb1,"model1_vaxgb.csv",row.names = F)

write.csv(vagam1,"model1_vagam.csv",row.names = F)

vaxgb1$variable[which(vaxgb1$variable=="AH1")] <- 'A'
vaxgb1$variable[which(vaxgb1$variable=="B")] <- 'C'
vaxgb1$variable[which(vaxgb1$variable=="o3")] <- 'E'
vaxgb1$variable[which(vaxgb1$variable=="no2")] <- 'F'
vaxgb1$variable[which(vaxgb1$variable=="pm10")] <- 'G'
vaxgb1$variable[which(vaxgb1$variable=="AH3")] <- 'B'

vagam1$variable[which(vagam1$variable=="AH1")] <- 'A'
vagam1$variable[which(vagam1$variable=="B")] <- 'C'
vagam1$variable[which(vagam1$variable=="o3")] <- 'E'
vagam1$variable[which(vagam1$variable=="no2")] <- 'F'
vagam1$variable[which(vagam1$variable=="pm10")] <- 'G'
vagam1$variable[which(vagam1$variable=="AH3")] <- 'B'



library(ggsci)

paxgb1 <- vaxgb1 %>% 
  group_by(year,age_group,variable) %>%
  summarise(yearly_attri=sum(death_attr)) %>%
  ggplot(data=)+
  geom_col(aes(x=as.factor(year),y=yearly_attri,fill=variable),position = "dodge")+
  labs(x='Year',y='Death counts')+
  theme(text=element_text(size=8))+
  facet_wrap(age_group ~ .,ncol = 2,scales = "free")+
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.key = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  scale_fill_manual(name = "Variables", 
                    labels = c('TypeA-H1','TypeA-H3','TypeB','O3','NO2','PM 10'),
                    values = c("#E64B35CC","#DF8F44CC","#F39B7FCC","#4DBBD5CC","#00A087CC","#79AF97CC"))




pagam1 <- vagam1 %>% 
  group_by(year,age_group,variable) %>%
  summarise(yearly_attri=sum(death_attr)) %>%
  ggplot(data=)+
  geom_col(aes(x=as.factor(year),y=yearly_attri,fill=variable),position = "dodge")+
  labs(x='Year',y='Death counts')+
  theme(text=element_text(size=8))+
  facet_wrap(age_group ~ .,ncol = 2,scales = "free")+
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.key = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  scale_fill_manual(name = "Variables", 
                    labels = c('TypeA-H1','TypeA-H3','TypeB','O3','NO2','PM 10'),
                    values = c("#E64B35CC","#DF8F44CC","#F39B7FCC","#4DBBD5CC","#00A087CC","#79AF97CC"))





#library(corrplot) # 载入corrplot包
library(Cairo)

#mcor <- cor(df[,-c(1,2,3)])
#colnames(mcor) <- c("Week","Temperture","Humidity","Year","TypeA_H1","TypeA_H3","TypeB","Adenovirus","RSV","O3","NO2","PM10")
#rownames(mcor) <- colnames(mcor)

#CairoPNG("corplot.png",width = 10000,height = 10000,dpi=1200)
#corrplot(mcor,type = {"upper"},tl.col = "black")
#dev.off()

CairoPNG("xgb_1.png",width = 10000,height = 10000,dpi=1200)
paxgb1
dev.off()

CairoPNG("gam_1.png",width = 10000,height = 10000,dpi=1200)
pagam1
dev.off()


plot_df <- plot_df1
trend1 <- ggplot(data=plot_df,aes(x=week_number))+
  geom_line(aes(y=actual_death,color='A'))+
  geom_line(data=subset(plot_df,week_number<=315),linetype=1,aes(y=predicted_xg,color='C'))+
  geom_line(data=subset(plot_df,week_number>315),linetype=6 ,aes(y=predicted_xg,color='C'))+
  geom_line(data=subset(plot_df,week_number<=315),linetype=1,aes(y=predicted_gam,color='B'))+
  geom_line(data=subset(plot_df,week_number>315),linetype=6 ,aes(y=predicted_gam,color='B'))+
  geom_vline(xintercept=315,linetype = 2)+
  scale_x_continuous(breaks = c(0, 54, 106, 158, 211, 263, 315, 367, 418),
                     labels = as.character(seq(from=2014,to=2022)))+
  labs(#title='Using 2010 to 2019 data to predict weekly all cause death counts of different age groups after 2020',
    x='Year',y='Death counts')+
  facet_grid(age_group ~ .,scales = "free")+
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.key = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  scale_color_manual(name = "Method", 
                     labels = c('Observed','GAM', 'XGBoost'),
                     values = c('grey20','blue2','red2'))+
  theme(text = element_text(size = 15) )

CairoPNG("model1.png",width = 6000,height = 3000,dpi=300)
trend1
dev.off()


time <- df %>% select(idx,year)
model1 <- merge(plot_df1,time,by.x = "week_number",by.y = "idx",all.x = T)
write.csv(model1,"result/model1.csv",row.names = F)

annual_excess <- model1 %>% 
  mutate(excess_gam=actual_death-predicted_gam,
         excess_xg=actual_death-predicted_xg) %>%
  group_by(year,age_group) %>% 
  summarise(xg=sum(excess_xg),gam=sum(excess_gam))

write.csv(annual_excess,'result/annual_excess.csv',quote=F,row.names = FALSE)


vaxgb1$variable[which(vaxgb1$variable=="A")] <- 'TypeA-H1'
vaxgb1$variable[which(vaxgb1$variable=="C")] <- 'TypeB'

vaxgb1$variable[which(vaxgb1$variable=="E")] <- 'O3'
vaxgb1$variable[which(vaxgb1$variable=="F")] <- 'NO2'
vaxgb1$variable[which(vaxgb1$variable=="G")] <- 'PM 10'
vaxgb1$variable[which(vaxgb1$variable=="B")] <- 'TypeA-H3'



vagam1$variable[which(vagam1$variable=="A")] <- 'TypeA-H1'
vagam1$variable[which(vagam1$variable=="C")] <- 'TypeB'

vagam1$variable[which(vagam1$variable=="E")] <- 'O3'
vagam1$variable[which(vagam1$variable=="F")] <- 'NO2'
vagam1$variable[which(vagam1$variable=="G")] <- 'PM 10'
vagam1$variable[which(vagam1$variable=="B")] <- 'TypeA-H3'



va_list <- c("vaxgb1","vagam1")

for (a in va_list){
  variable_attr <- get(a)
  
  bar_plot_df<-data.frame(category = character(),
                          excess_rate = numeric(),
                          upr = numeric(),
                          lwr = numeric(),
                          age_group = character(),
                          variable=character())
  variable_list<-unique(variable_attr$variable)
  age_list<-unique(variable_attr$age_group)
  for (i in variable_list){
    for (j in age_list){
      bar_plot_temp<-data.frame(category = c('Pre Pandemic','Pandemic'),
                                excess_rate = rep(NA,2),
                                upr = rep(NA,2),
                                lwr = rep(NA,2),
                                age_group = rep(j,2),
                                variable=rep(i,2))
      temp_df_pre <- variable_attr %>% filter(variable==i) %>% filter(age_group==j) %>%
        filter(year<2020)
      size <- dim(temp_df_pre)[1]
      excess_rate <- rep(NA,10000)
      for (n in 1:10000){
        temp_df<-temp_df_pre[sample(1:size,size,replace=TRUE),]
        excess_rate[n]<-sum(100*temp_df$death_attr/temp_df$pop)/6
      }
      
      bar_plot_temp$excess_rate[1]<-excess_rate[order(excess_rate)][5000]
      bar_plot_temp$lwr[1]<-excess_rate[order(excess_rate)][250]
      bar_plot_temp$upr[1]<-excess_rate[order(excess_rate)][9750]
      
      
      temp_df_pan<-variable_attr%>%filter(variable==i)%>%filter(age_group==j)%>%
        filter(year==2020 | year==2021)
      size<-dim(temp_df_pan)[1]
      excess_rate<-rep(NA,10000)
      for (n in 1:10000){
        temp_df<-temp_df_pan[sample(1:size,size,replace=TRUE),]
        excess_rate[n]<-sum(100*temp_df$death_attr/temp_df$pop)/2
      }
      
      bar_plot_temp$excess_rate[2]<-excess_rate[order(excess_rate)][5000]
      bar_plot_temp$lwr[2]<-excess_rate[order(excess_rate)][250]
      bar_plot_temp$upr[2]<-excess_rate[order(excess_rate)][9750]
      
      bar_plot_df<-rbind(bar_plot_df,bar_plot_temp)
    }
  }
  bar_plot_df$category<-factor(bar_plot_df$category,levels = c("Pre Pandemic","Pandemic"))
  
  paste <- paste("./result/box_",a,".png")
  CairoPNG(paste,width = 16000, height = 7000,dpi=1200)
  p <- ggplot(bar_plot_df,aes(x=variable,y=excess_rate,fill=category,group=category))+
    geom_bar(stat="identity",position = "dodge")+
    geom_errorbar(aes(ymin = lwr,ymax = upr),position=position_dodge())+
    facet_wrap(age_group~.,scales='free')+
    labs(x='Air pollutants and influenza proxies',y='Annual excess mortality rate (per 100,000 population)')+
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.key = element_blank(),
          strip.background = element_rect(colour="white", fill="white"))+
    theme(axis.text.x = element_text(angle = 45,hjust = 1) )
  print(p)
  dev.off()
  
}

#pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
# library("tidyverse")
# plot_df <- plot_df1
# wave_excess <- plot_df %>% 
#   mutate(wave = case_when(week_number >= 527 & week_number <= 534 ~ "1",
#                           week_number >= 535 & week_number <= 538 ~ "2",
#                           week_number >= 550 & week_number <= 561 ~ "3",
#                           week_number >= 569 & week_number <= 576 ~ "4")) %>% 
#            mutate(excess_gam=actual_death-predicted_gam,
#                   excess_xg=actual_death-predicted_xg) %>%
#            group_by(wave,age_group) %>% 
#            summarise(xg=sum(excess_xg),
#                      gam=sum(excess_gam))
#          write.csv(wave_excess,"result/wave_excess.csv",row.names = F)
#          
#          
#          wave_trend <- ggplot(data=subset(plot_df,week_number >= 524),aes(x=week_number))+
#            geom_line(aes(y=actual_death,color='Observed'))+
#            geom_line(aes(y=predicted_xg,color='XGBoost'))+
#            geom_line(aes(y=predicted_gam,color='GAM'))+
#            geom_vline(xintercept=527,linetype = 2)+
#            geom_vline(xintercept=534,linetype = 2)+
#            geom_vline(xintercept=538,linetype = 2)+
#            geom_vline(xintercept=550,linetype = 2)+
#            geom_vline(xintercept=561,linetype = 2)+
#            geom_vline(xintercept=569,linetype = 2)+
#            geom_vline(xintercept=576,linetype = 2)+
#            scale_x_continuous(breaks = c(527,534,538,550,561,569,576),
#                               labels = c("2020/01/28","2020/03/16","2020/04/12","2020/07/05","2020/09/21","2020/11/19","2021/01/08"))+
#            labs(#title='Using 2010 to 2019 data to predict weekly all cause death counts of different age groups after 2020',
#              x='Date',y='Death counts')+
#            facet_grid(age_group ~ .,scales = "free")+
#            theme_bw() +
#            theme(axis.line = element_line(color='black'),
#                  plot.background = element_blank(),
#                  panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank(),
#                  panel.border = element_blank(),
#                  legend.key = element_blank(),
#                  strip.background = element_rect(colour="white", fill="white"))+
#            scale_color_manual(name = "Method", 
#                               labels = c('Observed','GAM', 'XGBoost'),
#                               values = c('grey20','blue2','red2'))+
#            theme(axis.text.x = element_text(angle = 45,hjust = 1) )
#          
#          
#          
#          annual_excess<-data.frame(excess_death = numeric(),
#                                    upr = numeric(),
#                                    lwr = numeric(),
#                                    age_group = character(),
#                                    method=character())
#          age_list<-unique(plot_df$age_group)  
#          
#          for (i in age_list){
#            excess_temp<-data.frame(excess_death = rep(NA,2),
#                                    upr = rep(NA,2),
#                                    lwr = rep(NA,2),
#                                    age_group = rep(i,2),
#                                    method=c('gam','xg'))
#            df<-plot_df%>%filter(week_number >= 527 & week_number <= 534)%>%filter(age_group==i)%>%
#              mutate(excess=actual_death-predicted_gam)
#            excess_gam<-df$excess
#            ann<-rep(NA,10000)
#            for (j in 1:10000){
#              ann[j]<-sum(sample(excess_gam,length(excess_gam),replace=TRUE))
#            }
#            excess_temp$excess_death[1]<-ann[order(ann)][5000]
#            excess_temp$lwr[1]<-ann[order(ann)][250]
#            excess_temp$upr[1]<-ann[order(ann)][9750]
#            
#            
#            df<-plot_df%>%filter(week_number >= 527 & week_number <= 534)%>%filter(age_group==i)%>%
#              mutate(excess=actual_death-predicted_xg)
#            excess_xg<-df$excess
#            ann<-rep(NA,10000)
#            for (j in 1:10000){
#              ann[j]<-sum(sample(excess_xg,length(excess_xg),replace=TRUE))
#            }
#            excess_temp$excess_death[2]<-ann[order(ann)][5000]
#            excess_temp$lwr[2]<-ann[order(ann)][250]
#            excess_temp$upr[2]<-ann[order(ann)][9750]
#            annual_excess<-rbind(annual_excess,excess_temp)
#          }
#          
#          wave1 <- annual_excess %>% mutate(wave="1")
#          
#          annual_excess<-data.frame(excess_death = numeric(),
#                                    upr = numeric(),
#                                    lwr = numeric(),
#                                    age_group = character(),
#                                    method=character())
#          age_list<-unique(plot_df$age_group)  
#          
#          for (i in age_list){
#            excess_temp<-data.frame(excess_death = rep(NA,2),
#                                    upr = rep(NA,2),
#                                    lwr = rep(NA,2),
#                                    age_group = rep(i,2),
#                                    method=c('gam','xg'))
#            df<-plot_df%>%filter(week_number >= 535 & week_number <= 538)%>%filter(age_group==i)%>%
#              mutate(excess=actual_death-predicted_gam)
#            excess_gam<-df$excess
#            ann<-rep(NA,10000)
#            for (j in 1:10000){
#              ann[j]<-sum(sample(excess_gam,length(excess_gam),replace=TRUE))
#            }
#            excess_temp$excess_death[1]<-ann[order(ann)][5000]
#            excess_temp$lwr[1]<-ann[order(ann)][250]
#            excess_temp$upr[1]<-ann[order(ann)][9750]
#            
#            
#            df<-plot_df%>%filter(week_number >= 550 & week_number <= 561)%>%filter(age_group==i)%>%
#              mutate(excess=actual_death-predicted_xg)
#            excess_xg<-df$excess
#            ann<-rep(NA,10000)
#            for (j in 1:10000){
#              ann[j]<-sum(sample(excess_xg,length(excess_xg),replace=TRUE))
#            }
#            excess_temp$excess_death[2]<-ann[order(ann)][5000]
#            excess_temp$lwr[2]<-ann[order(ann)][250]
#            excess_temp$upr[2]<-ann[order(ann)][9750]
#            annual_excess<-rbind(annual_excess,excess_temp)
#          }
#          
#          wave2 <- annual_excess %>% mutate(wave="2")
#          
#          annual_excess<-data.frame(excess_death = numeric(),
#                                    upr = numeric(),
#                                    lwr = numeric(),
#                                    age_group = character(),
#                                    method=character())
#          age_list<-unique(plot_df$age_group)  
#          
#          for (i in age_list){
#            excess_temp<-data.frame(excess_death = rep(NA,2),
#                                    upr = rep(NA,2),
#                                    lwr = rep(NA,2),
#                                    age_group = rep(i,2),
#                                    method=c('gam','xg'))
#            df<-plot_df%>%filter(week_number >= 535 & week_number <= 538)%>%filter(age_group==i)%>%
#              mutate(excess=actual_death-predicted_gam)
#            excess_gam<-df$excess
#            ann<-rep(NA,10000)
#            for (j in 1:10000){
#              ann[j]<-sum(sample(excess_gam,length(excess_gam),replace=TRUE))
#            }
#            excess_temp$excess_death[1]<-ann[order(ann)][5000]
#            excess_temp$lwr[1]<-ann[order(ann)][250]
#            excess_temp$upr[1]<-ann[order(ann)][9750]
#            
#            
#            df<-plot_df%>%filter(week_number >= 569 & week_number <= 576)%>%filter(age_group==i)%>%
#              mutate(excess=actual_death-predicted_xg)
#            excess_xg<-df$excess
#            ann<-rep(NA,10000)
#            for (j in 1:10000){
#              ann[j]<-sum(sample(excess_xg,length(excess_xg),replace=TRUE))
#            }
#            excess_temp$excess_death[2]<-ann[order(ann)][5000]
#            excess_temp$lwr[2]<-ann[order(ann)][250]
#            excess_temp$upr[2]<-ann[order(ann)][9750]
#            annual_excess<-rbind(annual_excess,excess_temp)
#          }
#          
#          wave3 <- annual_excess %>% mutate(wave="3")
#          
#          annual_excess<-data.frame(excess_death = numeric(),
#                                    upr = numeric(),
#                                    lwr = numeric(),
#                                    age_group = character(),
#                                    method=character())
#          age_list<-unique(plot_df$age_group)  
#          
#          for (i in age_list){
#            excess_temp<-data.frame(excess_death = rep(NA,2),
#                                    upr = rep(NA,2),
#                                    lwr = rep(NA,2),
#                                    age_group = rep(i,2),
#                                    method=c('gam','xg'))
#            df<-plot_df%>%filter(week_number >= 535 & week_number <= 538)%>%filter(age_group==i)%>%
#              mutate(excess=actual_death-predicted_gam)
#            excess_gam<-df$excess
#            ann<-rep(NA,10000)
#            for (j in 1:10000){
#              ann[j]<-sum(sample(excess_gam,length(excess_gam),replace=TRUE))
#            }
#            excess_temp$excess_death[1]<-ann[order(ann)][5000]
#            excess_temp$lwr[1]<-ann[order(ann)][250]
#            excess_temp$upr[1]<-ann[order(ann)][9750]
#            
#            
#            df<-plot_df%>%filter(week_number >= 569 & week_number <= 576)%>%filter(age_group==i)%>%
#              mutate(excess=actual_death-predicted_xg)
#            excess_xg<-df$excess
#            ann<-rep(NA,10000)
#            for (j in 1:10000){
#              ann[j]<-sum(sample(excess_xg,length(excess_xg),replace=TRUE))
#            }
#            excess_temp$excess_death[2]<-ann[order(ann)][5000]
#            excess_temp$lwr[2]<-ann[order(ann)][250]
#            excess_temp$upr[2]<-ann[order(ann)][9750]
#            annual_excess<-rbind(annual_excess,excess_temp)
#          }
#          
#          wave4 <- annual_excess %>% mutate(wave="4")
#          
#          a <- rbind(wave1,wave2,wave3,wave4) %>% mutate(me="pi")
#          b <- cbind(b,a)
         
         