library(tidyverse)
library(Cairo)
library(corrplot)

load('all_cause_data/85+_weekly_mortality.rd')
load('all_cause_data/0-20_weekly_mortality.rd')
load('all_cause_data/20-40_weekly_mortality.rd')
load('all_cause_data/40-65_weekly_mortality.rd')
load('all_cause_data/65-85_weekly_mortality.rd')
load('all_cause_data/all_age_weekly_mortality.rd')
df_list<-c('20','40','65','85','max','all')
df_name<-c('0-19','20-39','40-64','65-84','85+','all')

temp_df<-data.frame(year = character(),
                    mortality_rate = numeric(),
                    age_group = character())



for (i in df_list){
  df<-get(paste('df',i,sep='_'))
  colnames(df)[2]<-'weekly_death'
  df$pop<-lowess(df$idx, df$pop)$y
  temp<-df%>%filter(year>2009)%>%mutate(mortality_rate=1000*weekly_death/pop,age_group=df_name[which(df_list==i)])%>%
    select(year,mortality_rate,age_group)
  temp$year[temp$year<2020]<-'Pre Pandemic'
  temp$year[temp$year==2020 | temp$year==2021]<-'Pandemic'
  temp_df<-rbind(temp_df,temp)
}

temp_df$disease<-'all-cause'



load('crd_data/85+_weekly_mortality_CRD.rd')
load('crd_data/0-20_weekly_mortality_CRD.rd')
load('crd_data/20-40_weekly_mortality_CRD.rd')
load('crd_data/40-65_weekly_mortality_CRD.rd')
load('crd_data/65-85_weekly_mortality_CRD.rd')
load('crd_data/all_age_weekly_mortality_CRD.rd')


temp_df_crd<-data.frame(year = character(),
                    mortality_rate = numeric(),
                    age_group = character())



for (i in df_list){
  df<-get(paste('df',i,sep='_'))
  colnames(df)[2]<-'weekly_death'
  df$pop<-lowess(df$idx, df$pop)$y
  temp<-df%>%filter(year>2009)%>%mutate(mortality_rate=1000*weekly_death/pop,age_group=df_name[which(df_list==i)])%>%
    select(year,mortality_rate,age_group)
  temp$year[temp$year<2020]<-'Pre Pandemic'
  temp$year[temp$year==2020 | temp$year==2021]<-'Pandemic'
  temp_df_crd<-rbind(temp_df_crd,temp)
}

temp_df_crd$disease<-'CRD'




load('PI_data/85+_weekly_mortality_PI.rd')
load('PI_data/0-20_weekly_mortality_PI.rd')
load('PI_data/20-40_weekly_mortality_PI.rd')
load('PI_data/40-65_weekly_mortality_PI.rd')
load('PI_data/65-85_weekly_mortality_PI.rd')
load('PI_data/all_age_weekly_mortality_PI.rd')


temp_df_PI<-data.frame(year = character(),
                        mortality_rate = numeric(),
                        age_group = character())



for (i in df_list){
  df<-get(paste('df',i,sep='_'))
  colnames(df)[2]<-'weekly_death'
  df$pop<-lowess(df$idx, df$pop)$y
  temp<-df%>%filter(year>2009)%>%mutate(mortality_rate=1000*weekly_death/pop,age_group=df_name[which(df_list==i)])%>%
    select(year,mortality_rate,age_group)
  temp$year[temp$year<2020]<-'Pre Pandemic'
  temp$year[temp$year==2020 | temp$year==2021]<-'Pandemic'
  temp_df_PI<-rbind(temp_df_PI,temp)
}

temp_df_PI$disease<-'PI'






temp<-rbind(temp_df,temp_df_crd,temp_df_PI)
colnames(temp)[1]<-'Category'
temp$Category<-factor(temp$Category,levels = c("Pre Pandemic","Pandemic"))
t <- ggplot(temp, aes(x=age_group, y=mortality_rate, fill=Category)) +
  geom_boxplot(notchwidth = 0.1)+
  scale_y_continuous(trans='log10')+
  facet_grid(. ~ disease,scales = "free")+
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.key = element_blank(),
        strip.background = element_rect(colour="white", fill="white"))+
  labs(x='Age Groups',y='Weekly mortality per million population')+
  theme(text = element_text(size = 15) )


#CairoPNG("boxplot.png",width = 16000, height = 7000,dpi=1200)
CairoPNG("boxplot.png",width = 5700,height = 3000,dpi=300)
t
dev.off()
