
# SPATIAL AND TEMPORAL VARIABILITY OF STEM EMISSIONS    (created by jbarba at 09072018)
## This script contents all data processing, analyses and plotting for stem emissions necessary for the hole study manuscript (except the direct calculations from SoilFluxPro)

# Libraries
  library(doBy)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(RColorBrewer)
  library(scales)
  library(nlme)
  library(lme4)
  library(svglite)
  library(dplyr)
  library(padr)
  library(zoo)
  library(lubridate)
  library(lmerTest)
  library(tidyverse)
  library(GGally)
  library(mgcv)
  library(mgcViz)
  library(gratia)
  library(visreg)
  library(sjPlot)
  library(merTools)
  library(ggResidpanel)
  library(emmeans)
  library(itsadug)
  library(cowplot)
  library(stargazer)

  setwd("C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/txt")

  
# 1) Data preparation ----
  
  # 1.1) Manual measurements ----
  
      # Merge LGR_field_notes with stem diameter, collar height, EVI, SF and GWT
      #If Atm Pressure was not merged with LGR_field_notes in the LGR_Flux_calculations.R for calculating the fluxes (and thus assumed constant Atm Press), Meteo data should be merged here.
      # LGR_fluxes.txt was created with the script LGR_Flux_calculations.R
      Manual_fluxes <- read.table("LGR_fluxes.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))    
      
      #Tree_diam: diameters for each tree at the heigh of each collar. Not necessary for calculating the fluxes
      Tree_diam <- read.table("SJ_tree_diameters.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
      #For_str: forest structure structure file for the plot. Not necessary for calculating the fluxes
      For_str <- read.table("SJ_forest_structure.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
      #GWT: ground water table measured every 15 min measured at the center of the plot. Not necessary for calculating the fluxes
      GWT<-read.table("SJ_GWL.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
      #SF: sap flow measured every 15 min at three different trees
      SF<-read.table("SJ_SF_full_experiment.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
      
      #  Manual_fluxes$Hour <- as.integer(sapply(strsplit(as.character(Manual_fluxes$Real_time), ":"), "[", 1))
      #  #Manual_fluxes (field notes) is a txt file with the information recorded in the field (LGR file name, code for samplig point, initial measurement time...[see the LGR_field_notes_example Excel])
      #  Manual_fluxes<-Manual_fluxes[complete.cases(Manual_fluxes), ]
      #  colnames(Manual_fluxes)[3]<-"Location"
      #  Manual_fluxes$DOY <- strftime(as.POSIXlt(Manual_fluxes$Date, format = "%m/%d/%Y"), format="%j")
      
      #This loop is for adding collar heights to Manual_fluxes
      for (i in 1:length (Manual_fluxes$Date)){
        for (k in 1:length (Tree_diam$Tree)){
          if (Manual_fluxes$Tree[i]==Tree_diam$Tree[k] & Manual_fluxes$Location[i]==Tree_diam$Location[k]) Manual_fluxes$Height[i]<-Tree_diam$Height[k]
        }
      }
      
      #This loop is for adding stem diameters to Manual_fluxes
      for (i in 1:length(Manual_fluxes$Date)){
        for (k in 1:length(For_str$Diam_est)){
          if (Manual_fluxes$Tree[i]==For_str$Tree_label[k]) {
            Manual_fluxes$Diam[i]<-For_str$Diam_est[k]
            Manual_fluxes$BA[i]<-For_str$BA[k]
          }  
        }
      }
      
      #This loop is for adding ground water table to Manual_fluxes
      GWT_mean <- summaryBy(Meters + Volts + Pulses ~ DOY + Hour + Month + Day + Date, data=GWT,FUN=c(mean,sd), na.rm=TRUE)
      
      for (i in 1:length(Manual_fluxes$Date)){
        for (k in 1:length(GWT_mean$Date)){
          if(Manual_fluxes$DOY[i]==GWT_mean$DOY[k] & Manual_fluxes$Hour[i]==GWT_mean$Hour[k]) Manual_fluxes$GWT[i]<-GWT_mean$Meters.mean[k]
        }
      }
      
      #This loop is for adding sap flow to Manual_fluxes
      SF_mean <- summaryBy(SF_avg ~ DOY + Hour, data=SF, FUN=c(mean,sd))
      SF_mean<-SF_mean[complete.cases(SF_mean), ]
      
      Manual_fluxes$SF<-NA
      for (i in 1:length(Manual_fluxes$Date)){
        for (k in 1:length(SF_mean$DOY)){
          if(Manual_fluxes$DOY[i]==SF_mean$DOY[k] & Manual_fluxes$Hour[i]==SF_mean$Hour[k]) Manual_fluxes$SF[i]<-SF_mean$SF_avg.mean[k] 
        }
      }
      
      ## This chunk of script was made by Andrew Hill from University of Delaware
      #This section is for creating Enhanced Vegetation Index (EVI) at dayly scale and adding it to Manual_fluxes
      
      #EVI Data
      #First go to: https://modis.ornl.gov/data.html register for a user name and password.  
      #Select global subsets tool, place marker on center of area of interest.
      #Select product, 16-day terra works well (MOD13Q1).
      #Sel Pixel size (0,0) will give min pixel size (250m x 250m).
      #Select date range.
      #View pixel footprint to insure good coverage of area of interest. 
      #Place order.
      #Recive email and access order.
      #Click CSV data. 
      #Choose product: statistics_250m_16_days_EVI.csv 
      #Open CSV and delete un-needed variables, leave Date and mean EVI.
      #Save CSV.
      
    #  EVI <- read.csv("EVI.csv", stringsAsFactors = FALSE)
   #   EVI$Date <- as_date(EVI$Date)
    #  EVI %>% pad(by="Date") -> EVI_open  #Open gaps
    #  EVI_open %>% mutate(EVI_interpol = na.approx(EVI_mean)) -> EVI_daily  #Interpolate open gaps at daily level.
      
   #   write.table(EVI_daily,file="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/txt/EVI_daily.txt",sep="\t")   
      
      #This loop is for EVI to Manual_fluxes
      EVI_daily$Date <- as_date(EVI_daily$Date)
      EVI_daily$DOY <- strftime(as.POSIXlt(EVI_daily$Date, format = "%y-%m-%d"), format="%j")
      
      for (i in 1:length(Manual_fluxes$Date)){
        for (k in 1:length(EVI_daily$Date)){
          if(Manual_fluxes$DOY[i]==EVI_daily$DOY[k]) Manual_fluxes$EVI[i]<-EVI_daily$EVI_interpol[k] 
        }
      }
      
      write.table(Manual_fluxes,file="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/txt/Manual_fluxes.txt",sep="\t")   

  # 1.2.) Automated measurements ----
      
      Stem_cont_flx<-read.table("SJ_continuous_stem_emissions.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
      date_time=strptime(Stem_cont_flx$Date,format='%m/%d/%Y')
      Stem_cont_flx$DOY<-(as.POSIXlt(date_time, format = "%d%b%y"))$yday
      Stem_cont_flx$DOYc<-Stem_cont_flx$DOY+as.numeric(Stem_cont_flx$Hour)/24+as.numeric(Stem_cont_flx$Minute)/(24*60)
      Stem_cont_flx$CH4_flux<-Stem_cont_flx$CH4_flux*1000  # Converts fluxes in nanomols
      Stem_cont_flx$N2O_flux<-Stem_cont_flx$N2O_flux*1000  # Converts fluxes in nanomols
      
      SWCv<-Stem_cont_flx$SWC^2*(-3.14E-7)+Stem_cont_flx$SWC*1.16E-3-0.612   #This equation is provided in the Licor manual for converting SWC from volts to volume units
      
      #This chunck is for assigning chambers to trees
      #Ports 1, 4 and 7 belong to soils
      #Ports 2, 5 and 8 belong to lower stems (+/- 50 cm)
      #Ports 3, 6 and 9 belong to upper stems (+/- 150 cm)
      for (i in 1:length(Stem_cont_flx$Date)){
        if(Stem_cont_flx$Port[i]==1|Stem_cont_flx$Port[i]==2|Stem_cont_flx$Port[i]==3) Stem_cont_flx$Tree[i]<-"1" 
        if(Stem_cont_flx$Port[i]==4|Stem_cont_flx$Port[i]==5|Stem_cont_flx$Port[i]==6) Stem_cont_flx$Tree[i]<-"2" 
        if(Stem_cont_flx$Port[i]==7|Stem_cont_flx$Port[i]==8|Stem_cont_flx$Port[i]==9) Stem_cont_flx$Tree[i]<-"3"
        if(Stem_cont_flx$Port[i]==1|Stem_cont_flx$Port[i]==4|Stem_cont_flx$Port[i]==7) Stem_cont_flx$Location[i]<-"1"
        if(Stem_cont_flx$Port[i]==2|Stem_cont_flx$Port[i]==5|Stem_cont_flx$Port[i]==8) Stem_cont_flx$Location[i]<-"2"
        if(Stem_cont_flx$Port[i]==3|Stem_cont_flx$Port[i]==6|Stem_cont_flx$Port[i]==9) Stem_cont_flx$Location[i]<-"3"
        if(Stem_cont_flx$Port[i]==1|Stem_cont_flx$Port[i]==4|Stem_cont_flx$Port[i]==7) Stem_cont_flx$code[i]<-"Soil"
        if(Stem_cont_flx$Port[i]==2|Stem_cont_flx$Port[i]==3|Stem_cont_flx$Port[i]==5|Stem_cont_flx$Port[i]==6|Stem_cont_flx$Port[i]==8|Stem_cont_flx$Port[i]==9) Stem_cont_flx$code[i]<-"Stem"
      }
      
      stem_h<-summaryBy(CO2_flux + CO2_R2 + CH4_flux + CH4_R2 + N2O_flux + N2O_R2 + Temperature + SWCv + DOYc~ DOY + Date + Year + Month + Day + Hour + Port + Location + Tree + code, data=Stem_cont_flx,FUN=c(mean), na.rm=TRUE)
      colnames(stem_h)[11:19]<-c("CO2_flux", "CO2_R2", "CH4_flux", "CH4_R2", "N2O_flux", "N2O_R2", "Temp", "SWC", "DOYc")
      
      #This loop is for adding collar heights to stem_h
      stem_h$Height<-NA
      for (i in 1:length(stem_h$Date)){
        for (k in 1:length (Tree_diam$Tree)){
          if (stem_h$Tree[i]==Tree_diam$Tree[k] & stem_h$Location[i]==as.character(Tree_diam$Location[k])) stem_h$Height[i]<-Tree_diam$Height[k]  
        }
      }
      
      #This loop is for adding stem diameters to stem_h
      stem_h$BA<-NA
      for (i in 1:length(stem_h$Date)){
        for (k in 1:length(For_str$Diam_est)){
          if (stem_h$Tree[i]==For_str$Tree_label[k]) {
            stem_h$Diam[i]<-For_str$Diam_est[k]
            stem_h$BA[i]<-For_str$BA[k]
          }  
        }
      }
      
      #This loop is for adding ground water table to stem_h
      GWT_mean <- summaryBy(Meters + Volts + Pulses ~ DOY + Hour + Month + Day + Date, data=GWT,FUN=c(mean,sd), na.rm=TRUE)
      
      stem_h$GWT<-NA
      for (i in 1:length(stem_h$Date)){
        for (k in 1:length(GWT_mean$Date)){
          if(stem_h$DOY[i]==GWT_mean$DOY[k] & stem_h$Hour[i]==GWT_mean$Hour[k]) stem_h$GWT[i]<-GWT_mean$Meters.mean[k]
        }
      }
      
      #This loop is for adding sap flow to stem_h
      SF_mean <- summaryBy(SF_avg ~ DOY + Hour, data=SF, FUN=c(mean,sd))
      SF_mean<-SF_mean[complete.cases(SF_mean), ]
      
      stem_h$SF<-NA
      for (i in 1:length(stem_h$Date)){
        for (k in 1:length(SF_mean$DOY)){
          if(stem_h$DOY[i]==SF_mean$DOY[k] & stem_h$Hour[i]==SF_mean$Hour[k]) stem_h$SF[i]<-SF_mean$SF_avg.mean[k] 
        }
      }
      
      write.table(stem_h,file="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/txt/SJ_continuous_stem_emissions_hourly.txt", sep="\t",row.names=F,col.names=T) 
      
      stem_d<-summaryBy(CO2_flux + CO2_R2 + CH4_flux + CH4_R2 + N2O_flux + N2O_R2 + Temp + SWC + Height + BA + Diam + GWT + SF ~ DOY + Date + Year + Month + Day + Port + Location + Tree + code, data=stem_h,FUN=c(mean,sd), na.rm=TRUE)
      colnames(stem_d)[10:35]<-c("CO2", "CO2_R2", "CH4", "CH4_R2", "N2O", "N2O_R2", "Temp", "SWC", "Height", "BA", "Diam", "GWT", "SF", "CO2sd", "CO2_R2sd", "CH4sd", "CH4_R2sd", "N2Osd", "N2O_R2sd", "Tempsd", "SWCsd", "Heightsd", "BAsd", "Diamsd", "GWTsd", "SFsd")
      
      write.table(stem_d,file="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/txt/SJ_continuous_stem_emissions_daily.txt", sep="\t",row.names=F,col.names=T) 
      
# 2) Figures ----
  
  # 2.1) Figure 1. Meteo + sap flow + water table ----
      
      #Soil temperature and SWC measured in the center of the plot, associated with the met station
      Meteo<-read.table("SJ_meteo_data.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
      
      F1_1<-ggplot(subset(Meteo, DOY<352&DOY>101), aes(x=DOYc)) +
        geom_point(aes(y=Soil_temp),na.rm=TRUE,size=0.7,color="red") +
        geom_point(aes(y=VWC*100-22),na.rm=TRUE, size=0.7) +
        scale_y_continuous(sec.axis = sec_axis(~ (. +22)/100, name=('SWC ('~m^3~m^-3*')'))) +
        ylab(bquote("Temperature  (" *degree*C~")")) +
        annotate("text", x=105, y=23, label= "a)",size=3.5, fontface =2) +
        #annotate("text", x=110, y=21, label= "SWC",size=4.5) +
        theme_bw() + theme(panel.grid = element_blank()) +
        theme(axis.text.x=element_blank(),
              axis.title.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title=element_text(size=10),
              axis.text.y = element_text(color = "red"),
              axis.title.y = element_text(color = "red"),
              axis.text.y.right = element_text(color = "black"),
              axis.title.y.right = element_text(color = "black"),
              plot.margin=unit(c(1,1,0,1), "cm"))
      
      #AtmPressure measured at 1.5m in the center of the plot
      F1_2<-ggplot(subset(Meteo, DOY>101&DOY<352), aes (x=DOYc, y=AtmPress)) +
        geom_line(na.rm=TRUE,size=0.7)+
        ylab(bquote("Atm Press (KPa)")) +
        ylim(98,104.5) +
        annotate("text", x=105, y=104, label= "b)",size=3.5, fontface =2) +
        theme_bw() + theme(panel.grid = element_blank()) +
        theme(axis.text.x=element_blank(),
              axis.title.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title=element_text(size=10),
              axis.text.y = element_text(color = "black"),
              axis.title.y = element_text(color = "black"),
              plot.margin=unit(c(-0.1,1,1,1), "cm"))
      
      
      #EVI (Enhanced vegetation index) remote-sensing measured in a patch of forest close to our plot. Data was measured every 14 days and extrapolated at daily basis)
      EVI_daily<-read.table("EVI_daily.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
      EVI_daily$DOY<-yday(EVI_daily$Date)
      
      F1_3<-ggplot(subset(EVI_daily, DOY>101&DOY<352), aes (x=DOY, y=EVI_interpol)) +
        geom_line(na.rm=TRUE,size=0.7)+
        ylab(bquote("EVI")) +
        annotate("text", x=105, y=0.72, label= "c)",size=3.5, fontface =2) +
        theme_bw() + theme(panel.grid = element_blank()) +
        theme(axis.text.x=element_blank(),
              axis.title.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y = element_text(color = "black"),
              axis.title.y = element_text(color = "black"),
              axis.title=element_text(size=10),
              plot.margin=unit(c(-1.1,1,1,1), "cm"))
      
      #Water table measured at the center of the plot
      F1_4<-ggplot(subset(GWT, DOY>101&DOY<352), aes (x=DOYc, y=Meters)) +
        geom_point(na.rm=TRUE,size=0.7)+
        xlab(bquote("DOY")) +
        ylab(bquote("GWT  (m)")) +
        annotate("text", x=105, y=-0.2, label= "d)",size=3.5, fontface =2) +
        geom_hline(yintercept = 0, linetype=2, size=0.5) +
        ylim(-1.5,0) +
        theme_bw() + theme(panel.grid = element_blank()) +
        theme(axis.text.x = element_text(color = "black"),
              axis.title.x = element_text(color = "black"),
              axis.text.y = element_text(color = "black"),
              axis.title.y = element_text(color = "black"),
              axis.title=element_text(size=10),
              plot.margin=unit(c(-1.1,1,1,1), "cm"))
      
      ## Combining all panels in of figure
      
      p1 <- ggplotGrob(F1_1)
      p2 <- ggplotGrob(F1_2)
      p3 <- ggplotGrob(F1_3)
      p4 <- ggplotGrob(F1_4)
      
      
      Meteo_graph<-grid.arrange(rbind(p1,p2,p3,p4, size = "first"))
      ggsave(filename="Meteo_graph.eps", plot=Meteo_graph, path="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/Figures", width=20, height=21.5)
      ggsave(filename="Meteo_graph.tiff", plot=Meteo_graph, path="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/Figures", width=20, height=21.5, units="cm", dpi=600)
      ggsave(filename="Meteo_graph.png", plot=Meteo_graph, path="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/Figures", width=15, height=17, units="cm", dpi=600)
      ##Saved as EPS, TIFF and PNG
      
  
  # 2.2) Figure 2. Manual stem CO2 emissions ----
  
   fluxes<-read.table("Manual_fluxes.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
  
    #CO2 flux for each tree and stem height with the same scale for all trees
    CO2_flux_tree_fixed<-ggplot(fluxes, aes(x=DOY, y=CO2_flux)) +
      geom_line(aes(color=as.factor(Location))) +
      geom_point(aes(color=as.factor(Location),shape=as.factor(Location))) +
      scale_color_brewer(palette="Dark2", labels = c("50cm","100cm","150cm"), name="Stem height") +  #Colors displayed for this palette: #1b9e77 for green, #d95f02 for orange and #7570b3 for purple
      scale_shape(labels=c("50cm","100cm","150cm"), name="Stem height") +
      geom_hline(yintercept = 0, linetype=2, size=0.5) +
      facet_wrap(~Tree) +
      xlab("DOY") + 
      ylab(bquote('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) +
      theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(colour = "grey20", size = 12),
            text = element_text(size = 16)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      theme(legend.position = c(1, 0),
            legend.justification = c(1.5, -0.2),
            legend.title = element_text(face =2),
            legend.text=element_text(size=12),
            axis.title=element_text(size=12),
            axis.text = element_text(color = "black", size=10)) +
      guides(color = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) 
    
    #CO2 flux for each tree and stem height with the adjusted scales for each tree
    CO2_flux_tree_free<-ggplot(fluxes, aes(x=DOY, y=CO2_flux)) +
      geom_line(aes(color=as.factor(Location))) +
      geom_point(aes(color=as.factor(Location),shape=as.factor(Location))) +
      scale_color_brewer(palette="Dark2", labels = c("50cm","100cm","150cm"), name="Stem height") +  #Colors displayed for this palette: #1b9e77 for green, #d95f02 for orange and #7570b3 for purple
      scale_shape(labels=c("50cm","100cm","150cm"), name="Stem height") +
      geom_hline(yintercept = 0, linetype=2, size=0.5) +
      facet_wrap(~Tree, scales="free") +
      xlab("DOY") + 
      ylab(bquote('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) +
      theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(colour = "grey20", size = 12),
            text = element_text(size = 16)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      theme(legend.position = c(1, 0),
            legend.justification = c(1.5, -0.2),
            legend.title = element_text(face =2),
            legend.text=element_text(size=12),
            axis.title=element_text(size=12),
            axis.text = element_text(color = "black", size=10)) +
      guides(color = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) 
    
    ggsave(file="CO2_manual_free.png", CO2_flux_tree_free, width=20, height=17, dpi=600, units='cm', device='png', path="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/Figures")
    ggsave(file="CO2_manual_free.tiff", CO2_flux_tree_free, width=20, height=17, dpi=600, units='cm', device='tiff', path="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/Figures")
    
    
  # 2.3) Figure 3. Manual stem CH4 emissions ----
  
    #CH4 flux for each tree and stem height with the same scale for all trees
    CH4_flux_tree_fixed<-ggplot(fluxes, aes(x=DOY, y=CH4_flux)) +
      geom_line(aes(color=as.factor(Location))) +
      geom_point(aes(color=as.factor(Location),shape=as.factor(Location))) +
      scale_color_brewer(palette="Dark2", labels = c("50cm","100cm","150cm"), name="Stem height") +  #Colors displayed for this palette: #1b9e77 for green, #d95f02 for orange and #7570b3 for purple
      scale_shape(labels=c("50cm","100cm","150cm"), name="Stem height") +
      geom_hline(yintercept = 0, linetype=2, size=0.5) +
      facet_wrap(~Tree) +
      xlab("DOY") + 
      ylab(bquote('CH'[4]~ 'flux  (nmol ' ~CH[4]~ m^-2~s^-1*')')) +
      theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(colour = "grey20", size = 12),
            text = element_text(size = 16)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      theme(legend.position = c(1, 0),
            legend.justification = c(1.5, -0.2),
            legend.title = element_text(face =2),
            legend.text=element_text(size=12),
            axis.title=element_text(size=12),
            axis.text = element_text(color = "black", size=10)) +
      guides(color = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) 
    
    #CH4 flux for each tree and stem height with the adjusted scales for each tree
    CH4_flux_tree_free<-ggplot(fluxes, aes(x=DOY, y=CH4_flux)) +
      geom_line(aes(color=as.factor(Location))) +
      geom_point(aes(color=as.factor(Location),shape=as.factor(Location))) +
      scale_color_brewer(palette="Dark2", labels = c("50cm","100cm","150cm"), name="Stem height") +  #Colors displayed for this palette: #1b9e77 for green, #d95f02 for orange and #7570b3 for purple
      scale_shape(labels=c("50cm","100cm","150cm"), name="Stem height") +
      geom_hline(yintercept = 0, linetype=2, size=0.5) +
      facet_wrap(~Tree, scales="free") +
      xlab("DOY") + 
      ylab(bquote('CH'[4]~ 'flux  (nmol ' ~CH[4]~ m^-2~s^-1*')')) +
      theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(colour = "grey20", size = 12),
            text = element_text(size = 16)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      theme(legend.position = c(1, 0),
            legend.justification = c(1.5, -0.2),
            legend.title = element_text(face =2),
            legend.text=element_text(size=12),
            axis.title=element_text(size=12),
            axis.text = element_text(color = "black", size=10)) +
      guides(color = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) 
  
    ggsave(file="CH4_manual_free.png", CH4_flux_tree_free, width=20, height=17, dpi=600, units='cm', device='png', path="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/Figures")
    ggsave(file="CH4_manual_free.tiff", CH4_flux_tree_free, width=20, height=17, dpi=600, units='cm', device='tiff', path="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/Figures")
    
 
  # 2.4)  Figure 4. Variance partitioning ----
    # The aim of this section is split the variance of the whole manual data between trees, location within trees, residual (temporal variability included within residual).
    
    mod.var.CO2 <-lmer(CO2_flux_t ~ 1 + (1|Tree/Location) + (1|DOY), data=fluxes)
    summary(mod.var.CO2)
    
    mod.var.CH4 <-lmer(CH4_flux_bc_z ~ 1 + (1|Tree/Location) + (1|DOY), data=fluxes)
    summary(mod.var.CH4)
    
    split_var<-data.frame()
    split_var[1,1]<-"Location"
    split_var[2,1]<-"DOY"
    split_var[3,1]<-"Tree"
    split_var[4,1]<-"Residual"
    split_var[5,1]<-"Location"
    split_var[6,1]<-"DOY"
    split_var[7,1]<-"Tree"
    split_var[8,1]<-"Residual"
    split_var[1:4,2]<-"CO2"
    split_var[5:8,2]<-"CH4"
    split_var[1,3]<-as.data.frame(VarCorr(mod.var.CO2))[1,'vcov']
    split_var[2,3]<-as.data.frame(VarCorr(mod.var.CO2))[2,'vcov']
    split_var[3,3]<-as.data.frame(VarCorr(mod.var.CO2))[3,'vcov']
    split_var[4,3]<-as.data.frame(VarCorr(mod.var.CO2))[4,'vcov']
    split_var[5,3]<-as.data.frame(VarCorr(mod.var.CH4))[1,'vcov']
    split_var[6,3]<-as.data.frame(VarCorr(mod.var.CH4))[2,'vcov']
    split_var[7,3]<-as.data.frame(VarCorr(mod.var.CH4))[3,'vcov']
    split_var[8,3]<-as.data.frame(VarCorr(mod.var.CH4))[4,'vcov']
    split_var[1,4]<-as.data.frame(VarCorr(mod.var.CO2))[1,'sdcor']
    split_var[2,4]<-as.data.frame(VarCorr(mod.var.CO2))[2,'sdcor']
    split_var[3,4]<-as.data.frame(VarCorr(mod.var.CO2))[3,'sdcor']
    split_var[4,4]<-as.data.frame(VarCorr(mod.var.CO2))[4,'sdcor']
    split_var[5,4]<-as.data.frame(VarCorr(mod.var.CH4))[1,'sdcor']
    split_var[6,4]<-as.data.frame(VarCorr(mod.var.CH4))[2,'sdcor']
    split_var[7,4]<-as.data.frame(VarCorr(mod.var.CH4))[3,'sdcor']
    split_var[8,4]<-as.data.frame(VarCorr(mod.var.CH4))[4,'sdcor']
    names(split_var)<-c("Component","Gas","Variance","sd")
    
    split_var$Gas <- factor(split_var$Gas, ordered = TRUE, levels = c("CO2", "CH4"))
    split_var$Component <- factor(split_var$Component, ordered = TRUE, levels = c("Residual", "DOY", "Location", "Tree"))
    
    x_axis <- c(expression(""~CO[2]),
                expression(""~CH[4]))
    
    plot_variance<-ggplot(split_var, aes(x = Gas, y=Variance*100, fill=Component)) +
      geom_bar(stat = 'identity') +
      scale_fill_brewer(palette="Dark2") +
      ylab("Variance (%)") +
      scale_x_discrete(labels= x_axis) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            legend.text=element_text(size=10),
            axis.title=element_text(size=10),
            axis.text = element_text(color = "black", size=10)) +
      guides(fill = guide_legend(reverse = TRUE))
    
    ggsave(file="Variance_partitioning.png", plot_variance, width=10, height=10, dpi=600, units='cm', device='png', path="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/Figures")
    ggsave(file="Variance_partitioning.tiff", plot_variance, width=10, height=10, dpi=600, units='cm', device='tiff', path="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/Figures")
    
    
    
  # 2.5) Figure 5. Automated stem emissions (CO2 & CH4) ----
    
    CO2_cont_1<-ggplot(subset(stem_d,code=="Stem"&Tree=="1"), aes(DOY, CO2, col=as.factor(Port), shape=as.factor(Port))) +
      geom_ribbon(aes(ymin = CO2 - CO2sd, ymax = CO2 + CO2sd, fill=as.factor(Port)), alpha=0.4, colour=NA) +
      scale_fill_manual(values=c("#1b9e77","#7570b3")) +
      geom_point(na.rm=TRUE, size=0.7) +
      geom_line(na.rm=TRUE, size=0.5) +
      scale_color_manual(values=c("#1b9e77","#7570b3")) + 
      geom_hline(yintercept = 0, linetype=2, size=0.5) +
      ylim(-2,12.5) +
      annotate("text", x = 285, y = 12, label = sprintf('\u25B2'), col="#7570b3", size=3) +
      annotate("text", x = 285, y = 10.5, label = sprintf('\u25CF'), col="#1b9e77", size=4) +
      ylab(bquote('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) +
      annotate("text", x=320, y=12, label= "150 cm",color="#7570b3",size=4, fontface =2) +
      annotate("text", x=320, y=10.5, label= "50 cm",color="#1b9e77",size=4, fontface =2) +
      annotate("text", x=103, y=12, label= "a)",size=4.5, fontface =2) +
      annotate("text", x=225, y=12, label= "Tree 1",size=4.5, fontface =2) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=12),
            axis.text = element_text(color = "black", size=12))
    
    CO2_cont_2<-ggplot(subset(stem_d,code=="Stem"&Tree=="2"), aes(DOY, CO2, col=as.factor(Port), shape=as.factor(Port))) +
      geom_ribbon(aes(ymin = CO2 - CO2sd, ymax = CO2 + CO2sd, fill=as.factor(Port)), alpha=0.4, colour=NA) +
      scale_fill_manual(values=c("#1b9e77","#7570b3")) +
      geom_point(na.rm=TRUE, size=0.7) +
      geom_line(na.rm=TRUE, size=0.5) +
      scale_color_manual(values=c("#1b9e77","#7570b3")) +  
      geom_hline(yintercept = 0, linetype=2, size=0.5) +
      ylim(-2,12.5) +
      ylab(bquote('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) +
      annotate("text", x=103, y=12, label= "c)",size=4.5, fontface =2) +
      annotate("text", x=225, y=12, label= "Tree 2",size=4.5, fontface =2) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=12),
            axis.text = element_text(color = "black", size=12))
    
    CO2_cont_3<-ggplot(subset(stem_d,code=="Stem"&Tree=="3"), aes(DOY, CO2, col=as.factor(Port), shape=as.factor(Port))) +
      geom_ribbon(aes(ymin = CO2 - CO2sd, ymax = CO2 + CO2sd, fill=as.factor(Port)), alpha=0.4, colour=NA) +
      scale_fill_manual(values=c("#1b9e77","#7570b3")) +
      geom_point(na.rm=TRUE, size=0.7) +
      geom_line(na.rm=TRUE, size=0.5) +
      scale_color_manual(values=c("#1b9e77","#7570b3")) +  
      geom_hline(yintercept = 0, linetype=2, size=0.5) +
      ylim(-2,12.5) +
      ylab(bquote('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) +
      annotate("text", x=103, y=12, label= "e)",size=4.5, fontface =2) +
      annotate("text", x=225, y=12, label= "Tree 3",size=4.5, fontface =2) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=12),
            axis.text = element_text(color = "black", size=12)) 
    
    #Upper CH4 fluxes for tree 1 are so hight that it is hard to see seasonal patterns for Lower stem if both are plotted together
    #If we want to plot both chambers in different scales or using secondary axes, we will need the next chunk of script.
    for (i in 1:length(stem_d$Date)){
      if(stem_d$Port[i]==2) {
        stem_d$Lower_tree_1[i]<-stem_d$CH4[i] 
        stem_d$Lower_tree_1sd[i]<-stem_d$CH4sd[i]
      }
      else{
        stem_d$Lower_tree_1[i]<-NA 
        stem_d$Lower_tree_1sd[i]<-NA
      }
      if(stem_d$Port[i]==3) {
        stem_d$Upper_tree_1[i]<-stem_d$CH4[i]
        stem_d$Upper_tree_1sd[i]<-stem_d$CH4sd[i]
      }
      else{
        stem_d$Upper_tree_1[i]<-NA 
        stem_d$Upper_tree_1sd[i]<-NA
      }
      if(stem_d$Port[i]!=3) {
        stem_d$Other[i]<-stem_d$CH4[i]
        stem_d$Othersd[i]<-stem_d$CH4sd[i]
      }
      else{
        stem_d$Other[i]<-NA 
        stem_d$Othersd[i]<-NA
      }
    }
    
    CH4_cont_1<-ggplot(subset(stem_d,code=="Stem"&Tree=="1"), aes(DOY, CH4, col=as.factor(Port), shape=as.factor(Port))) +
      geom_ribbon(aes(ymin = CH4 - CH4sd, ymax = CH4 + CH4sd, fill=as.factor(Port)), alpha=0.4, colour=NA) +
      scale_fill_manual(values=c("#1b9e77","#7570b3")) +
      geom_point(na.rm=TRUE, size=0.7) +
      geom_line(na.rm=TRUE, size=0.5) +
      scale_color_manual(values=c("#1b9e77","#7570b3")) + 
      geom_hline(yintercept = 0, linetype=2, size=0.5) +
      ylab(bquote('CH'[4]~ 'flux  (nmol ' ~CH[4]~ m^-2~s^-1*')')) +
      annotate("text", x=103, y=36, label= "b)",size=4.5, fontface =2) +
      annotate("text", x=225, y=36, label= "Tree 1",size=4.5, fontface =2) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=12),
            axis.text = element_text(color = "black", size=12))    
    
    CH4_cont_2<-ggplot(subset(stem_d,code=="Stem"&Tree=="2"), aes(DOY, CH4, col=as.factor(Port), shape=as.factor(Port))) +
      geom_ribbon(aes(ymin = CH4 - CH4sd, ymax = CH4 + CH4sd, fill=as.factor(Port)), alpha=0.4, colour=NA) +
      scale_fill_manual(values=c("#1b9e77","#7570b3")) +
      geom_point(na.rm=TRUE, size=0.7) +
      geom_line(na.rm=TRUE, size=0.5) +
      scale_color_manual(values=c("#1b9e77","#7570b3")) + 
      geom_hline(yintercept = 0, linetype=2, size=0.5) +
      ylim(-2,4) +
      #xlab("Date") + 
      ylab(bquote('CH'[4]~ 'flux  (nmol ' ~CH[4]~ m^-2~s^-1*')')) +
      annotate("text", x=103, y=3.8, label= "d)",size=4.5, fontface =2) +
      annotate("text", x=225, y=3.8, label= "Tree 2",size=4.5, fontface =2) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            #axis.text.x=element_blank(),
            #axis.title.x=element_blank(),
            #axis.ticks.x=element_blank(),
            axis.title=element_text(size=12),
            axis.text = element_text(color = "black", size=12))    
    
    CH4_cont_3<-ggplot(subset(stem_d,code=="Stem"&Tree=="3"), aes(DOY, CH4, col=as.factor(Port), shape=as.factor(Port))) +
      geom_ribbon(aes(ymin = CH4 - CH4sd, ymax = CH4 + CH4sd, fill=as.factor(Port)), alpha=0.4, colour=NA) +
      scale_fill_manual(values=c("#1b9e77","#7570b3")) +
      geom_point(na.rm=TRUE, size=0.7) +
      geom_line(na.rm=TRUE, size=0.5) +
      scale_color_manual(values=c("#1b9e77","#7570b3")) +  
      geom_hline(yintercept = 0, linetype=2, size=0.5) +
      ylim(-2,4) +
      ylab(bquote('CH'[4]~ 'flux  (nmol ' ~CH[4]~ m^-2~s^-1*')')) +
      annotate("text", x=103, y=3.8, label= "f)",size=4.5, fontface =2) +
      annotate("text", x=225, y=3.8, label= "Tree 3",size=4.5, fontface =2) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=12),
            axis.text = element_text(color = "black", size=12))   
    
    CO2_CH4_cont <- grid.arrange(CO2_cont_1, CH4_cont_1, CO2_cont_2, CH4_cont_2, CO2_cont_3, CH4_cont_3, ncol=2)
    
    ggsave(file="GHG_cont_time_series.png", CO2_CH4_cont, width=20, height=20, dpi=600, units='cm', device='png', path="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/Figures")
    ggsave(file="GHG_cont_time_series.tiff", CO2_CH4_cont, width=20, height=20, dpi=600, units='cm', device='tiff', path="C:/Users/barbafej/Dropbox/SJ_stem_fluxes/SJ_studies/SJ_full_paper/Figures")
  
      
  # 2.6) Figure 6. Density plot and mean fluxes of CO2 and CH4 ----
        
        ggplot(subset(stem_h,code=="Stem"),aes(CH4_flux,DOYc))+geom_point()
    
    density_CH4<-ggplot(data=subset(fluxes, CH4_flux<20), aes(x=CH4_flux,..scaled..)) + ## Do I have to keep the ..scaled..?
          geom_density(adjust=1,colour=F,fill="red",alpha=0.5)+
          geom_density(data=subset(stem_h,code=="Stem"),aes(x=CH4_flux,..scaled..),adjust=1,colour=F,fill="dodgerblue4",alpha=0.5)+
      #ggplot(A,aes(x,y)) +geom_point() +geom_point(data=B,colour='red') + xlim(0, 10)    
      #geom_density(aes(CH4c,..scaled..),adjust=2,fill="red",colour=F,alpha=0.5)+
          #scale_y_continuous(expand=c(0,0), limits=c(0,1.1))+  
          xlim(-10,25)+
          ylab(bquote("Density")) +
          xlab(bquote('CH'[4]~ 'flux  (nmol ' ~CH[4]~ m^-2~s^-1*')')) +
          #geom_segment(aes(x = mean_c, y = 0, xend = mean_c, yend = 0.67),linetype=2, size=0.5,color="black") +
          #geom_segment(aes(x = mean_m, y = 0, xend = mean_m, yend = 0.93),linetype=2, size=0.5,color="black") +
          theme_bw() + theme(panel.grid = element_blank()) +
          theme(axis.text.x=element_text(color="black",size=12),
                axis.text.y=element_text(color="black",size=12),
                axis.title=element_text(size=12),
                axis.ticks.x=element_line(color="black"),
                axis.ticks.y=element_line(color="black"))
          
    density_CO2<-ggplot(fluxes, aes(x=CO2_flux,..scaled..)) + ## Do I have to keep the ..scaled..?
      geom_density(adjust=1,colour=F,fill="red",alpha=0.5)+
      geom_density(data=subset(stem_h,code=="Stem"),aes(x=CO2_flux,..scaled..),adjust=1,colour=F,fill="dodgerblue4",alpha=0.5)+
      #ggplot(A,aes(x,y)) +geom_point() +geom_point(data=B,colour='red') + xlim(0, 10)    
      #geom_density(aes(CH4c,..scaled..),adjust=2,fill="red",colour=F,alpha=0.5)+
      #scale_y_continuous(expand=c(0,0), limits=c(0,1.1))+  
      xlim(-10,25)+
      ylab(bquote("Density")) +
      xlab(bquote('CH'[4]~ 'flux  (nmol ' ~CH[4]~ m^-2~s^-1*')')) +
      #geom_segment(aes(x = mean_c, y = 0, xend = mean_c, yend = 0.67),linetype=2, size=0.5,color="black") +
      #geom_segment(aes(x = mean_m, y = 0, xend = mean_m, yend = 0.93),linetype=2, size=0.5,color="black") +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(axis.text.x=element_text(color="black",size=12),
            axis.text.y=element_text(color="black",size=12),
            axis.title=element_text(size=12),
            axis.ticks.x=element_line(color="black"),
            axis.ticks.y=element_line(color="black"))     
    
  # 2.7) Figure 7. Stem and soil CO2 and CH4 concentrations ----
      
    # 2.7 a & b Stem concentrations ----  
  
        ##ggplot2, doBy grifExtra are required for running this chunk
        Stem_c<-read.table("SJ_stem_concentrations.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
        Stem_c$Location <- factor(Stem_c$Location, ordered = TRUE, levels = c("1", "2", "3"))    
        Stem_c$CH4t<-log10(Stem_c$CH4)
        Stem_c$CO2t<-log10(Stem_c$CO2)
        Stem_c$code<-paste0(Stem_c$Tree,Stem_c$Location)
        Stem_c[46,2]<-"1"
        Stem_c[48,2]<-"1"
        Stem_c[50,2]<-"1"
        
        #This is an ANOVA with repeated measures (Tree as an error factor) for testing differences between Heights  
        modelAOV <- aov(CH4t~Location+Error(factor(Tree)), data = Stem_c)
        summary(modelAOV)
        model_CO2_lm <- nlme::lme(CO2t~Location,random=~1|Tree,data=Stem_c)
        summary(model_CO2_lm)
        anova(model_CO2_lm)
        plot(model_CO2_lm)
        
        model_CH4_lm <- nlme::lme(CH4t~Location,random=~1|Tree,data=Stem_c)
        summary(model_CH4_lm)
        anova(model_CH4_lm)
        plot(model_CH4_lm)
        
        modelAOV <- aov(CO2t~Location+Error(factor(Tree)), data = Stem_c)
        summary(modelAOV)
        
        ## Plot [CH4] and [CO2] at 08/20/2018
        
        Stem_c$std_height<-as.numeric(Stem_c$Location)*50
        
        Stem_c$Location<-as.factor(Stem_c$Location)
        
        CO2_axis <- c("10^3", "10^4", "10^5")
        levels(CO2_axis) <- c("10^3", "10^4", "10^5")
        
        CO2_stem_concentration<-ggplot(Stem_c, aes(x=std_height, y=CO2t, fill=factor(std_height))) +
          geom_violin(trim=F) +
          geom_point() +
          scale_fill_brewer(palette="Dark2") +
          scale_x_continuous(limits=c(5,180), breaks=c(0, 50, 100, 150)) + 
          scale_y_continuous(position="right", limits=c(2.5,5.5), breaks=c(3,4,5), labels=parse(text = levels(CO2_axis))) + 
          coord_flip() +
          annotate("text", x=15, y=5.3, label= "a)",size=4.5, fontface =2) +
          xlab('Stem Height (cm)') + 
          ylab('[CO'[2]~']  (ppmv) (log scale)') +
          theme_bw() + theme(panel.grid = element_blank()) +
          theme(legend.position="none",
                axis.title=element_text(size=12),
                axis.text = element_text(color = "black", size=12),
                plot.margin=unit(c(0.2,0,0.2,0.2), "cm"))
        
        CH4_axis <- c("10[-2]","1", "10[2]", "10[4]", "10[6]")
        levels(CH4_axis) <- c("10^-2","1", "10^2", "10^4", "10^6")
        
        CH4_stem_concentration<-ggplot(Stem_c, aes(x=std_height, y=CH4t, fill=factor(std_height))) +
          geom_violin(trim=F) +
          geom_point() +
          scale_fill_brewer(palette="Dark2") +  #Colors displayed for this palette: #1b9e77 for green, #d95f02 for orange and #7570b3 for purple
          scale_x_continuous(limits=c(5,180), breaks=c(0, 50, 100, 150)) + 
          scale_y_continuous(position="right", breaks=c(-2,0,2,4,6), labels=parse(text = levels(CH4_axis))) + 
          coord_flip() +
          annotate("text", x=15, y=6.8, label= "b)",size=4.5, fontface =2) +
          xlab('Stem Height (cm)')  +
          ylab('[CH'[4]~'] (ppmv) (log scale)') +
          theme_bw() + theme(panel.grid = element_blank()) +
          theme(legend.position="none",
                axis.title.x=element_text(size=12),
                axis.title.y=element_blank(),
                axis.text.x = element_text(color = "black", size=12),
                axis.text.y = element_blank(),
                axis.ticks.y=element_blank(),
                plot.margin=unit(c(0.2,0.2,0.2,-0.1), "cm"))
        
        stem_concentrations <- grid.arrange(CO2_stem_concentration,CH4_stem_concentration, nrow=1) 
        ggsave(file="Stem_concentrations.png", stem_concentrations, width=16, height=10, dpi=600, units='cm', device='png', path="//anr.udel.edu/files/shares/jbarba/My Documents/SJ_Data/SJ_stem_fluxes/SJ_full_paper/Figures")
        ggsave(file="Stem_concentrations.tiff", stem_concentrations, width=16, height=10, dpi=600, units='cm', device='tiff', path="//anr.udel.edu/files/shares/jbarba/My Documents/SJ_Data/SJ_stem_fluxes/SJ_full_paper/Figures")
        ggsave(file="Stem_concentrations.svg", plot=stem_concentrations, width=16, height=10, units='cm', path="//anr.udel.edu/files/shares/jbarba/My Documents/SJ_Data/SJ_stem_fluxes/SJ_full_paper/Figures")
        
    # 2.7 c & d Soil concentrations ----
        
        Soil_concentrations<-read.table("SJ_soil_concentrations.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
        Soil_concentrations<-na.omit(Soil_concentrations)
        
        Soil_c<-summaryBy(CH4 + CO2 ~ Date + Depth, data=Soil_concentrations,FUN=c(mean,sd), na.rm=TRUE)
         
        ## Plot [CH4] and [CO2] measured on 10/06/2017, 08/20/2018 and 02/31/2019
        
        CO2_soil<-ggplot(data=Soil_c, aes(x=Depth, y=CO2.mean, color=as.factor(Date))) +
          geom_point(data=Soil_concentrations, aes(x=Depth, y=CO2, color=as.factor(Date)), size=1.5) +
          geom_line(na.rm=TRUE, size=0.75) +
          scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a")) +
          scale_x_reverse(limits = c(160, -3),breaks = c(0, 10, 25, 50, 75, 100, 150), label = c("0", "10", "25", "50", "75", "100", "GW")) +
          scale_y_continuous(position = "top", limits=c(0,12500)) +
          coord_flip() +
          geom_hline(yintercept = Soil_c[8,4], linetype=2, size=0.5) +
          xlab('Soil depth (cm)') + 
          ylab('[CO'[2]~']  (ppmv)') +
          annotate("text", x=0, y=11750, label= "c)",size=4.5, fontface =2) +
          annotate("text", x=120, y=10000, label= "October 2017",size=4, fontface =2, color="#4daf4a") +
          annotate("text", x=135, y=10000, label= "August 2018",size=4, fontface =2, color="#377eb8") +
          annotate("text", x=150, y=10000, label= "March 2019",size=4, fontface =2, color="#e41a1c") +
          theme_bw() + theme(panel.grid = element_blank()) +
          theme(legend.position="none",
                axis.title.x=element_text(size=12),
                axis.title.y=element_blank(),
                axis.text.x = element_text(color = "black", size=12),
                axis.text.y = element_blank(),
                axis.ticks.y=element_blank(),
                plot.margin=unit(c(0.2,0,0.2,0.2), "cm"))
        
        CH4_soil<-ggplot(data=Soil_c, aes(x=Depth, y=CH4.mean, color=as.factor(Date))) +
          geom_point(data=Soil_concentrations, aes(x=Depth, y=CH4, color=as.factor(Date)), size=1.5) +
          geom_line(na.rm=TRUE, size=0.75) +
          scale_color_manual(values=c("#e41a1c","#377eb8","#4daf4a")) +
          scale_x_reverse(limits = c(160, -3),breaks = c(0, 10, 25, 50, 75, 100, 150), label = c("0", "10", "25", "50", "75", "100", "GW"), position="left") +
          scale_y_continuous(position = "top", limits=c(-0.1,5.3)) +
          coord_flip() +
          geom_hline(yintercept = Soil_c[8,3], linetype=2, size=0.5) +
          xlab('Soil depth (cm)') + 
          ylab('[CH'[4]~']  (ppmv)') +
          annotate("text", x=0, y=5, label= "d)",size=4.5, fontface =2) +
          theme_bw() + theme(panel.grid = element_blank()) +
          theme(legend.position="none",
                axis.title.x=element_text(size=12),
                axis.title.y=element_blank(),
                axis.text.x = element_text(color = "black", size=12),
                axis.text.y = element_blank(),
                plot.margin=unit(c(0.2,0.2,0.2,-0.1), "cm")) 
        
        soil_concentrations<-grid.arrange(CO2_soil, CH4_soil, ncol=2)
        ggsave(file="Soil_concentrations.png", soil_concentrations, width=16, height=10, dpi=600, units='cm', device='png', path="//anr.udel.edu/files/shares/jbarba/My Documents/SJ_Data/SJ_stem_fluxes/SJ_full_paper/Figures")
        ggsave(file="Soil_concentrations.tiff", soil_concentrations, width=16, height=10, dpi=600, units='cm', device='tiff', path="//anr.udel.edu/files/shares/jbarba/My Documents/SJ_Data/SJ_stem_fluxes/SJ_full_paper/Figures")
        ggsave(file="Soil_concentrations.svg", plot=soil_concentrations, width=16, height=10, units='cm', path="//anr.udel.edu/files/shares/jbarba/My Documents/SJ_Data/SJ_stem_fluxes/SJ_full_paper/Figures")
 
  
  # 2.8) Tree cores incubations ----
  
        incubations<-read.table("SJ_cores_incubations.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
        
        # Shapiro-Wilk normality test for Heartwood
        with(incubations, shapiro.test(log(CO2[Tissue == "Heartwood"])))
        # Shapiro-Wilk normality test for Sapwood
        with(incubations, shapiro.test(log(CO2[Tissue == "Sapwood"]))) 
        #Shapiro test indicates that Heartwood and Sapwood data are normally distributed (p>0.05)
        var.test(log(CO2) ~ Tissue, data = incubations)
        #The p-value of F-test is p < 0.001. It's smaller than the significance level alpha = 0.05, so there are differences between the variances of the two sets of data. 
        #Therefore, we cannot use the classic t-test witch assume equality of the two variances.
        
        t.test(CO2 ~ Tissue, data=incubations, alternative="greater", var.equal = FALSE)
        
        
        
        # Shapiro-Wilk normality test for Heartwood
        with(incubations, shapiro.test(CH4[Tissue == "Heartwood"]))
        # Shapiro-Wilk normality test for Sapwood
        with(incubations, shapiro.test(CH4[Tissue == "Sapwood"])) 
        #Shapiro test indicates that Heartwood and Sapwood data are normally distributed (p>0.05)
        var.test(CH4 ~ Tissue, data = incubations)
        #The p-value of F-test is p = 0.1171. It's greater than the significance level alpha = 0.05. 
        #In conclusion, there is no significant difference between the variances of the two sets of data. 
        #Therefore, we can use the classic t-test witch assume equality of the two variances.
        
        t.test(CH4 ~ Tissue, data=incubations, alternative="greater", var.equal = FALSE)
        
        incubations$Tree<- stringr::str_sub(incubations$Sample,1,1)
        incubations$Location<- stringr::str_sub(incubations$Sample,2,2) 
        
        model.ch4.1<- nlme::lme(CH4 ~ Tissue*Location, data=incubations,random=~1|Tree)
        summary(model.ch4.1)
        plot(model.ch4.1)
        anova(model.ch4.1)
        
        model.ch4.2<- nlme::lme(CH4 ~ Tissue, data=incubations,random=~1|Tree/Location)
        summary(model.ch4.2)
        plot(model.ch4.2)
        anova(model.ch4.2)
        
        model.ch4.3<- nlme::lme(CH4 ~ Tree*Tissue, data=incubations,random=~1|Location)
        summary(model.ch4.3)
        plot(model.ch4.3)
        anova(model.ch4.3)
                
        inc_CO2<-ggplot(incubations, aes(x=Tissue, y=CO2, fill=factor(Tissue))) +
          geom_violin(trim=F) +
          scale_fill_manual(values=c("#fc8d59","#91cf60")) +
          geom_point(size=0.9) +
          ylab('CO'[2]~' (ppmv)') +
          annotate("text", x=2.35, y=1075, label= "a)",size=4.5, fontface =2) +
          annotate("text", x=1.5, y=975, label= "***",size=6, fontface =2) +
          theme_bw() + theme(panel.grid = element_blank()) +
          theme(legend.position="none",
                axis.text.x=element_blank(),
                axis.title=element_text(size=13),
                axis.title.x=element_blank(),
                axis.text.y = element_text(color = "black", size=12),
                plot.margin=unit(c(0.2,0.2,0.2,0.2), "cm")) 
        
        inc_CH4<-ggplot(incubations, aes(x=Tissue, y=CH4, fill=factor(Tissue))) +
          geom_violin(trim=F) +
          scale_fill_manual(values=c("#fc8d59","#91cf60")) +
          geom_point(size=0.9) +
          ylab('CH'[4]~' (ppmv)') +
          scale_y_continuous(limits=c(180,330), breaks=c(175,225,275,325)) + 
          annotate("text", x=2.35, y=325, label= "b)",size=4.5, fontface =2) +
          annotate("text", x=1.5, y=315, label= "***",size=6, fontface =2) +
          theme_bw() + theme(panel.grid = element_blank()) +
          theme(legend.position="none",
                axis.text.x=element_text(color = "black", size=13),
                axis.title=element_text(size=13),
                axis.title.x=element_blank(),
                axis.text.y = element_text(color = "black", size=12),
                plot.margin=unit(c(-0.3,0.2,0.2,0.2), "cm"))
        
        inc <- grid.arrange(inc_CO2,inc_CH4, ncol=1) 
        ggsave(file="Core_incubations.png", inc, width=8, height=16, dpi=600, units='cm', device='png', path="//anr.udel.edu/files/shares/jbarba/My Documents/SJ_Data/SJ_stem_fluxes/SJ_full_paper/Figures")
        ggsave(file="Core_incubations.tiff", inc, width=8, height=16, dpi=600, units='cm', device='tiff', path="//anr.udel.edu/files/shares/jbarba/My Documents/SJ_Data/SJ_stem_fluxes/SJ_full_paper/Figures")


# 3) Data analysis -----    
    fluxes<-read.table("Manual_fluxes.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
    
    #Data preparation    
    fluxes$Location <- factor(fluxes$Location, ordered = TRUE, levels = c("1", "2", "3"))
    fluxes$Locationf <- factor(fluxes$Location)
    fluxes<-subset(fluxes,CH4_flux<10) #There is one CH4 manual measurement with a very high value, which is claerly an statistical outlier. 
    # Box-cox transformation of CH4
    fluxes$CH4_flux_bc <- car::bcnPower(fluxes$CH4_flux,
                                lambda=car::powerTransform(fluxes$CH4_flux,family='bcnPower')[['lambda']],
                                gamma=car::powerTransform(fluxes$CH4_flux,family='bcnPower')[['gamma']])  
    fluxes$CH4_flux_bc_z <- scale(fluxes$CH4_flux_bc,scale=TRUE)
    fluxes$Temp_z <- scale(fluxes$Temp,scale=TRUE)
    fluxes$VWC_z <- scale(fluxes$VWC,scale=TRUE)
    fluxes$Diam_z <- scale(fluxes$Diam,scale=TRUE)
    fluxes$Temp_c <- scale(fluxes$Temp,scale=FALSE)
    fluxes$VWC_c <- scale(fluxes$VWC,scale=FALSE)
    fluxes$Diam_c <- scale(fluxes$Diam,scale=FALSE)
    
    fluxes.cont<-read.table("SJ_continuous_stem_emissions_daily.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
    fluxes.cont$Port<-as.factor(fluxes.cont$Port)
    fluxes.cont$Tree<-as.factor(fluxes.cont$Tree)
    fluxes.cont<-subset(fluxes.cont, fluxes.cont$code=="Stem")
    fluxes.cont$Locationf <-factor(fluxes.cont$Location)
    fluxes.cont$Location <- factor(fluxes.cont$Location, ordered = TRUE, levels = c("2", "3"))
      
    
    fluxes.cont$CH4_bc <- car::bcnPower(fluxes.cont$CH4,
                           lambda=car::powerTransform(fluxes.cont$CH4,family='bcnPower')[['lambda']],
                           gamma=car::powerTransform(fluxes.cont$CH4,family='bcnPower')[['gamma']])
    fluxes.cont$CH4_bc_z = scale(fluxes.cont$CH4_bc,scale=TRUE)
    fluxes.cont$Temp_z = scale(fluxes.cont$Temp,scale=TRUE)
    fluxes.cont$SWC_z = scale(fluxes.cont$SWC,scale=TRUE)
    fluxes.cont$Temp_c = scale(fluxes.cont$Temp,scale=FALSE)
    fluxes.cont$SWC_c = scale(fluxes.cont$SWC,scale=FALSE)
    
  #3.1) Manual measurements of CO2 stem emissions ----       
    
    # model without temporal autocorrelation
    lmm_co2_lme<- lme(log(CO2_flux)~Temp*Location+VWC*Location+Diam,
                      random=~Temp|Tree,data=fluxes)
    # model with temporal autocorrelation
    lmm_co2_AR1 <- update(lmm_co2_lme,cor=corAR1(form=~1|Tree/Location))
    
    anova(lmm_co2_lme,lmm_co2_AR1)
    
    par(mfcol=c(1,2))
    acf(residuals(lmm_co2_lme,type='normalized'))
    acf(residuals(lmm_co2_AR1,type='normalized'))
    
    plot(lmm_co2_AR1)
    summary(lmm_co2_AR1)
    MuMIn::r.squaredGLMM(lmm_co2_AR1)
    car::Anova(lmm_co2_AR1,3)
    
  #3.2) Manual measurements of CH4 stem emissions ----  
    
    # untransformed ch4 flux
    lmm_ch4_lme<- lme(CH4_flux~Temp*Location+VWC*Location+Diam,
                      random=~1|Tree,data=fluxes)
    plot(lmm_ch4_lme) #residuals don't llok good
    
    # use box-cox transformed ch4 fluxt
    lmm_ch4bc_lme<- lme(CH4_flux_bc~Temp*Location+VWC*Location+Diam,
                        random=~1|Tree,data=fluxes)
    plot(lmm_ch4bc_lme)
    
    # add corAR1
    lmm_ch4bc_AR1 <- update(lmm_ch4bc_lme,
                            cor=corAR1(form=~DOY|Tree/Location))
    plot(lmm_ch4bc_AR1)
    anova(lmm_ch4bc_lme,lmm_ch4bc_AR1)
    
    par(mfcol=c(1,2))
    acf(residuals(lmm_ch4bc_lme,type='normalized'))
    acf(residuals(lmm_ch4bc_AR1,type='normalized'))
    
    plot(lmm_ch4bc_AR1)
    summary(lmm_ch4bc_AR1)
    MuMIn::r.squaredGLMM(lmm_ch4bc_AR1)
    car::Anova(lmm_ch4bc_AR1,3)
    
  #3.3) Automated CO2 measurements ----
    
    # model without temporal autocorrelation
    co2.lme <- lme(log(CO2)~Temp*Locationf*SWC, random= ~Temp|Tree,
                   na.action=na.omit,data=fluxes.cont)
    
    # model with temporal autocorrelation
    co2.lme.AR1 <- update(co2.lme,cor=corAR1(form=~DOY|Tree/Locationf))
    
    par(mfcol=c(1,2))
    acf(residuals(co2.lme,type='normalized'))
    acf(residuals(co2.lme.AR1,type='normalized'))
    
    plot(co2.lme.AR1)
    summary(co2.lme.AR1)
    MuMIn::r.squaredGLMM(co2.lme.AR1)
    car::Anova(co2.lme.AR1,3)
    
    cowplot::plot_grid(nrow=2,
                       sjPlot::plot_model(co2.lme.AR1,type='pred',terms=c('SWC','Locationf')),
                       sjPlot::plot_model(co2.lme.AR1,type='pred',terms=c('Temp','Locationf'))
    )
    
  #3.3) Automated CH4 measurements ----
    # model without temporal autocorrelation
    ch4bc.lme<- lme(CH4_bc~Temp*Locationf*SWC, random= ~Temp|Tree,
                    data=fluxes.cont,na.action=na.omit)
    
    # model variance as a function of Tree to consider heteroscedasticity
    ch4bc.lme.varIdent<- lme(CH4_bc~Temp*Locationf*SWC, random= ~Temp|Tree,
                             weights=varIdent(form=~1|Tree),
                             data=fluxes.cont,na.action=na.omit)
    # Add corAR1
    ch4bc.lme.varIdentAR1<- update(ch4bc.lme.varIdent,cor=corAR1(form=~DOY|Tree/Locationf))
    
    par(mfcol=c(1,2))
    acf(residuals(ch4bc.lme.varIdent,type='normalized'))
    acf(residuals(ch4bc.lme.varIdentAR1,type='normalized'))
    
    anova(ch4bc.lme.varIdent,ch4bc.lme.varIdentAR1)
    
    plot(ch4bc.lme.varIdentAR1)
    summary(ch4bc.lme.varIdentAR1)
    MuMIn::r.squaredGLMM(ch4bc.lme.varIdentAR1)
    car::Anova(ch4bc.lme.varIdentAR1,3)
    
    sjPlot::plot_model(ch4bc.lme.varIdentAR1,type='pred',terms=c('SWC','Locationf'))
    
    
    
    
    
    # model without temporal autocorrelation
    co2.lme <- lme(log(CO2)~Temp*Location*SWC, random= ~Temp|Tree,
                   na.action=na.omit,data=fluxes.cont)
    
    # model with temporal autocorrelation
    co2.lme.AR1 <- update(co2.lme,cor=corAR1(form=~DOY|Tree/Location))
    
    anova(co2.lme,co2.lme.AR1)
    
    par(mfcol=c(1,2))
    acf(residuals(co2.lme,type='normalized'))
    acf(residuals(co2.lme.AR1,type='normalized'))
    
    plot(co2.lme.AR1)
    summary(co2.lme.AR1)
    MuMIn::r.squaredGLMM(co2.lme.AR1)
    car::Anova(co2.lme.AR1,3)
    
  #3.4) Automated CH4 measurements ----
    
    # model without temporal autocorrelation
    ch4bc.lme<- lme(CH4_bc~Temp*Location*SWC, random= ~Temp|Tree,
                    data=fluxes.cont,na.action=na.omit)
    
    # model variance as a function of Tree to consider heteroscedasticity
    ch4bc.lme.varIdent<- lme(CH4_bc~Temp*Location*SWC, random= ~Temp|Tree,
                             weights=varIdent(form=~1|Tree),
                             data=fluxes.cont,na.action=na.omit)
    # Add corAR1
    ch4bc.lme.varIdentAR1<- update(ch4bc.lme.varIdent,cor=corAR1(form=~DOY|Tree/Location))
    
    anova (ch4bc.lme.varIdent,ch4bc.lme.varIdentAR1)
    par(mfcol=c(1,2))
    acf(residuals(ch4bc.lme.varIdent,type='normalized'))
    acf(residuals(ch4bc.lme.varIdentAR1,type='normalized'))
    
    plot(ch4bc.lme.varIdentAR1)
    summary(ch4bc.lme.varIdentAR1)
    MuMIn::r.squaredGLMM(ch4bc.lme.varIdentAR1)
    anova(ch4bc.lme.varIdentAR1)
    
          
  #3.1 Relation between CO2 and CH4 fluxes. There are large differences in magnitudes between trees and locations (different heights), so we put   
      # This is the model you suggested
      
      mod.CH4_CO2<-lmer(CH4_flux ~ CO2_flux * Location + (Location|Tree), data=fluxes)
      
      #Warning message:
      #  In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
      #                 Model failed to converge with max|grad| = 0.00638437 (tol = 0.002, component 1)
      
      fluxes$CH4_flux_t<-scale(fluxes$CH4_flux)
      fluxes$CO2_flux_t<-scale(fluxes$CO2_flux)
      mod.CH4_CO2t<-lmer(CH4_flux_t ~ CO2_flux_t * Location + (Location|Tree), data=fluxes)
      #boundary (singular) fit: see ?isSingular
      
      summary(mod.CH4_CO2t)
      #I don't know how to interpret this summary. There is no p-value for the whole model of for the variables, as it used to be in the nlme models  
    
      # Additionally, if I'm interested in the relation between CH4 and CO2, why we are not using the // CH4_flux_t ~ CO2_flux_t + (Location|Tree) //? 
    
  #3.2 Relation between CO2 or CH4 fluxes with stem height
      fluxes<-read.table("Manual_fluxes.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
      fluxes$Location <- factor(fluxes$Location, ordered = TRUE, levels = c("1", "2", "3"))
      fluxes<-subset(fluxes,CH4_flux<10) #There is one CH4 manual measurement with a very high value, which is claerly an statistical outlier. 
    #CH4 with Location (relative position with stem height (factor))  
      #This is an ANOVA with repeated measures (Tree as an error factor) for testing differences between Heights  
      # We need the library(lmerTest) to obtain the p-values
    library(multcomp)
      #En Rafa diu que el segon model te molt mes sentit
      # mod.CO2_Location<-lmer(CO2_flux ~ Location + (1|Tree/Location) + (1|DOY), data=fluxes) 

      mod.CO2_loc_diam<-lmer(log(CO2_flux) ~ Location + Diam + (1|Tree) + (1|DOY), data=fluxes)
      summary(mod.CO2_loc_diam)
      anova(mod.CO2_loc_diam)
      confint(mod.CO2_loc_diam)
      summary(glht(mod.CO2_loc_diam, linfct=mcp(Location="Tukey")), test = adjusted(type = "bonferroni"))
    fluxes$Location <- factor(fluxes$Location, ordered=FALSE)     
      mod.CH4_loc_diam<-lmer(CH4_flux_bc_z ~ Location + Diam + (1|Tree) + (1|DOY), data=fluxes)
      summary(mod.CH4_loc_diam)
      anova(mod.CH4_loc_diam)
      confint(mod.CH4_loc_diam)
      confint.merMod(mod.CH4_loc_diam,oldNames=FALSE)
      summary(glht(mod.CH4_loc_diam, linfct=mcp(Location="Tukey")), test = adjusted(type = "bonferroni"))
      
      modelAOV <- aov(CH4_flux~Location + Error(factor(Tree)), data = fluxes)
      summary(modelAOV)
      TukeyHSD(modelAOV, "Location", ordered = TRUE) #Error in UseMethod("TukeyHSD") : no applicable method for 'TukeyHSD' applied to an object of class "c('aovlist', 'listof')"
      plot(TukeyHSD(modelAOV, "Location")) #Error in UseMethod("TukeyHSD") : no applicable method for 'TukeyHSD' applied to an object of class "c('aovlist', 'listof')"
      
      
      
      mod.CO2_diam<-lmer(CO2_flux ~ Diam + (1|Tree) + (1|DOY), data=fluxes)
      summary(mod.CO2_diam)
      anova(mod.CO2_diam)

      
      mod.CH4_diam<-lmer(CH4_flux ~ Diam + (1|Tree) + (1|DOY), data=fluxes)
      summary(mod.CH4_diam)
      anova(mod.CH4_diam)
      
      
      
      
      #I cannot run the TukeyHSD so I run the next lme model
      mod.CH4_loc <- lme(CH4_flux ~ Location, random = ~ 1 | Tree, 
                     weights=varIdent(~1|Location), #You put this line of code months ago, but I don't have any idea what does it mean (or do)
                     data=fluxes)
      summary(mod.CH4_loc)
      plot(mod.CH4_loc)
      plot(fluxes$Location, residuals(mod.CH4_Location, type = "norm"))
      anova(mod.CH4_loc)
      
      require(multcomp)
      summary(glht(mod.CH4_Location, linfct=mcp(Location="Tukey")), test = adjusted(type = "bonferroni"))
      
    #CO2 with Location (relative position with stem height (factor))  
      #This is an ANOVA with repeated measures (Tree as an error factor) for testing differences between Heights  
      
      modelAOV <- aov(CO2_flux~Location + Error(factor(Tree)), data = fluxes)
      summary(modelAOV)
      TukeyHSD(modelAOV, "Location", ordered = TRUE) #Error in UseMethod("TukeyHSD") : no applicable method for 'TukeyHSD' applied to an object of class "c('aovlist', 'listof')"
      plot(TukeyHSD(modelAOV, "Location")) #Error in UseMethod("TukeyHSD") : no applicable method for 'TukeyHSD' applied to an object of class "c('aovlist', 'listof')"
      
      #I cannot run the TukeyHSD so I run the next lme model
      mod.CO2_loc <- lme(CO2_flux ~ Location, random = ~ 1 | Tree, 
                         weights=varIdent(~1|Location), #You put this line of code months ago, but I don't have any idea what does it mean (or do)
                         data=fluxes)
      summary(mod.CO2_loc)
      plot(mod.CO2_loc) #This plot doesn't llok so nice
      plot(fluxes$Location, residuals(mod.CO2_loc, type = "norm"))
      anova(mod.CO2_loc)
      
      summary(glht(mod.CH4_loc, linfct=mcp(Location="Tukey")), test = adjusted(type = "bonferroni"))
      
  
# model time series ch4. I have to put all the other variables and rename the models
    
    modch4.randint<- nlme::lme(CH4_flux~(Location+Temp)^2,data=fluxes,
                               na.action=na.omit,
                               random=~1|Tree,
                               correlation=nlme::corAR1(),
                               weights=varIdent(~1|Location)) #You put this line of code months ago, but I don't have any idea what does it mean (or do)
                               
    
    summary(modch4.randint)
    anova(modch4.randint)
    visreg::visreg(modch4.randint,'Temp','Location')
    
    
    modco2.randint<- nlme::lme(CO2_flux~(Location+Temp)^2,data=fluxes,
                               na.action=na.omit,
                               random=~1|Tree,
                               correlation=nlme::corAR1(),
                               weights=varIdent(~1|Location))
    
    summary(modco2.randint)
    anova(modco2.randint)
    visreg::visreg(modco2.randint,'Temp','Location')
    

## Mean fluxes
    #Mean, sd and 95% CI for automated measurements
      stems<-subset(Stem_cont_flx, Port==2 | Port==3 | Port==5 | Port==6 | Port==8 | Port==9)  
      stems_by_tree<-summaryBy(CO2_flux + CH4_flux ~ Tree, data=stems,FUN=mean, na.rm=TRUE)
      mean_CO2_aut<-mean(stems_by_tree$CO2_flux, na.rm=TRUE)
      sd_CO2_aut<-sd(stems_by_tree$CO2_flux, na.rm=TRUE)
      CI_CO2_aut<-qnorm(0.975)*sd_CO2_aut/sqrt(length(stems_by_tree$CO2_flux))
      mean_CH4_aut<-mean(stems_by_tree$CH4_flux*1000, na.rm=TRUE)
      sd_CH4_aut<-sd(stems_by_tree$CH4_flux*1000, na.rm=TRUE)
      CI_CH4_aut<-qnorm(0.975)*sd_CH4_aut/sqrt(length(stems_by_tree$CH4_flux))
    #Mean manual measurements
      fluxes_by_tree<-summaryBy(CO2_flux + CH4_flux ~ Tree, data=fluxes,FUN=mean, na.rm=TRUE)
      mean_CO2_manual<-mean(fluxes_by_tree$CO2_flux, na.rm=TRUE)
      sd_CO2_manual<-sd(fluxes_by_tree$CO2_flux, na.rm=TRUE)
      CI_CO2_manual<-qnorm(0.975)*sd_CO2_manual/sqrt(length(fluxes_by_tree$CO2_flux))
      mean_CH4_manual<-mean(fluxes_by_tree$CH4_flux, na.rm=TRUE)
      sd_CH4_manual<-sd(fluxes_by_tree$CH4_flux, na.rm=TRUE)
      CI_CH4_manual<-qnorm(0.975)*sd_CH4_manual/sqrt(length(fluxes_by_tree$CH4_flux))
    #% of measurements with net uptake
      uptake_automated<-length(subset(stems,CH4_flux<0)$CH4_flux)/length(stems$CH4_flux[!is.na(stems$CH4_flux)])*100
      uptake_manual<-length(subset(fluxes,CH4_flux<0)$CH4_flux)/length(fluxes$CH4_flux[!is.na(fluxes$CH4_flux)])*100
    #Mean, sd and 95% CI for manual measurements BY STEM HEIGHT
      fluxes_by_location<-summaryBy(CO2_flux + CH4_flux ~ Tree + Location, data=fluxes,FUN=mean, na.rm=TRUE)
      fluxes_by_location2<-summaryBy(CO2_flux.mean + CH4_flux.mean ~ Location, data=fluxes_by_location,FUN=c(mean,sd), na.rm=TRUE)
      
      
     #Mean, sd and 95% CI for automated measurements BY STEM HEIGHT
      stems_by_location<-summaryBy(CO2_flux + CH4_flux ~ Tree + Location, data=stems,FUN=mean, na.rm=TRUE)
      stems_by_location2<-summaryBy(CO2_flux + CH4_flux ~ Location, data=stems,FUN=c(mean,sd), na.rm=TRUE)
      
      #mean_CO2_aut_location<-mean(stems_by_location$CO2_flux, na.rm=TRUE)
      #sd_CO2_aut_location<-sd(stems_by_location$CO2_flux, na.rm=TRUE)
      CI_CO2_aut_50<-qnorm(0.975)*stems_by_location2[1,4]/sqrt(length(unique(stems_by_location$Tree)))
      CI_CO2_aut_150<-qnorm(0.975)*stems_by_location2[2,4]/sqrt(length(unique(stems_by_location$Tree)))
      
      mean_CH4_aut_location<-mean(stems_by_location$CH4_flux*1000, na.rm=TRUE)
      sd_CH4_aut_location<-sd(stems_by_location$CH4_flux*1000, na.rm=TRUE)
      CI_CH4_aut_location<-qnorm(0.975)*sd_CH4_aut_location/sqrt(length(stems_by_location$CH4_flux))
          
    #CO2 and CH4 vs. stem heigth

      
  # 3.5) CO2 - CH4 relationship ----
   
      #lmm_co2chc_lme<- lme(CO2_flux ~ CH4_flux,
      #                    random=~1|Tree,data=fluxes)
      #Continuous data
      co2_ch4.lme <- lme(CH4_bc~CO2*Locationf, random= ~1|Tree/Locationf,
                     na.action=na.omit,cor=corAR1(form=~DOY|Tree/Locationf),data=fluxes.cont)
      
      
      CO2_CH4_automated<-ggplot(fluxes.cont, aes(x=CO2, y=CH4_bc)) +
        #geom_line(aes(color=as.factor(Location))) +
        geom_point() +
        #scale_color_brewer(palette="Dark2", labels = c("50cm","100cm","150cm"), name="Stem height") +  #Colors displayed for this palette: #1b9e77 for green, #d95f02 for orange and #7570b3 for purple
        #scale_shape(labels=c("50cm","100cm","150cm"), name="Stem height") +
        #geom_hline(yintercept = 0, linetype=2, size=0.5) +
        #facet_wrap(~Tree) +
        #xlab("DOY") + 
        xlab(bquote('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) +
        ylab(bquote('CH'[4]~ 'flux  (Box-Cox transformed)')) +
        theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
              axis.text.y = element_text(colour = "grey20", size = 12),
              text = element_text(size = 16)) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        theme(axis.title=element_text(size=12),
              axis.text = element_text(color = "black", size=10))  
      
      # model with temporal autocorrelation
      #co2_ch4.lme.AR1 <- update(co2_ch4.lme,cor=corAR1(form=~DOY|Tree/Locationf))
      
      par(mfcol=c(1,1))
      acf(residuals(co2_ch4.lme,type='normalized'))
      acf(residuals(co2_ch4.lme,type='normalized'))
      
      plot(co2_ch4.lme)
      summary(co2_ch4.lme)
      MuMIn::r.squaredGLMM(co2_ch4.lme)
      car::Anova(co2_ch4.lme,3)
      
      sjPlot::plot_model(co2_ch4.lme,type='pred',terms=c('CO2','Locationf'))
      
      #Manual data
      CO2_CH4_manual<- lme(CH4_flux_bc~log(CO2_flux)*Locationf,
                        random=~1|Tree/Locationf,cor=corAR1(form=~1|Tree/Locationf),data=fluxes)
      # model with temporal autocorrelation
      lmm_co2_AR1 <- update(CO2_CH4_manual,cor=corAR1(form=~1|Tree/Location))
      
      anova(lmm_co2_lme,lmm_co2_AR1)
      
      par(mfcol=c(1,2))
      acf(residuals(CO2_CH4_manual,type='normalized'))
      acf(residuals(lmm_co2_AR1,type='normalized'))
      
      plot(CO2_CH4_manual)
      summary(CO2_CH4_manual)
      MuMIn::r.squaredGLMM(CO2_CH4_manual)
      car::Anova(CO2_CH4_manual,3)
      
      sjPlot::plot_model(CO2_CH4_manual,type='pred',terms=c('CO2_flux','Locationf'))
      
      CO2_CH4<-ggplot(fluxes, aes(x=log(CO2_flux), y=CH4_flux_bc)) +
        #geom_line(aes(color=as.factor(Location))) +
        geom_point() +
        #scale_color_brewer(palette="Dark2", labels = c("50cm","100cm","150cm"), name="Stem height") +  #Colors displayed for this palette: #1b9e77 for green, #d95f02 for orange and #7570b3 for purple
        #scale_shape(labels=c("50cm","100cm","150cm"), name="Stem height") +
        #geom_hline(yintercept = 0, linetype=2, size=0.5) +
        #facet_wrap(~Tree) +
        #xlab("DOY") + 
        #ylab(bquote('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')) +
        theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
              axis.text.y = element_text(colour = "grey20", size = 12),
              text = element_text(size = 16)) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        theme(legend.position = c(1, 0),
              legend.justification = c(1.5, -0.2),
              legend.title = element_text(face =2),
              legend.text=element_text(size=12),
              axis.title=element_text(size=12),
              axis.text = element_text(color = "black", size=10)) +
        guides(color = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) 
      
# SUPPLEMENTARY MATERIAL ----
    # We have sapwood and heartwood density, length and moisture for tree cores sampled at 150cm. We are re-running the GLMs for manual measurements (subset 150cm fluxes) to see is there is any effect of thoise variables
    #These analyses will be included in the suplementary material
    
      wood<-read.table("Wood_properties.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
      
      #fluxes_WP is a subset of fluxes for measurements at 150 cm with the wood properties (density, length and moisture)
      fluxes_wp<-subset(fluxes,Location=="3" & Tree!="J") #I don't have wood data for tree J

      fluxes_wp$CH4_flux_bc <- car::bcnPower(fluxes_wp$CH4_flux,
                                          lambda=car::powerTransform(fluxes_wp$CH4_flux,family='bcnPower')[['lambda']],
                                          gamma=car::powerTransform(fluxes_wp$CH4_flux,family='bcnPower')[['gamma']])  
      fluxes_wp$CH4_flux_bc_z <- scale(fluxes_wp$CH4_flux_bc,scale=TRUE)
      fluxes_wp$Temp_z <- scale(fluxes_wp$Temp,scale=TRUE)
      fluxes_wp$VWC_z <- scale(fluxes_wp$VWC,scale=TRUE)
      fluxes_wp$Diam_z <- scale(fluxes_wp$Diam,scale=TRUE)
      fluxes_wp$Temp_c <- scale(fluxes_wp$Temp,scale=FALSE)
      fluxes_wp$VWC_c <- scale(fluxes_wp$VWC,scale=FALSE)
      fluxes_wp$Diam_c <- scale(fluxes_wp$Diam,scale=FALSE)
      
      for (t in 1:length(fluxes_wp$Date)){
        for (g in 1:length(wood$Tree_Name)){
          if (fluxes_wp$Tree[t]==wood$Tree_Name[g]) {
            fluxes_wp$Density[t]<- wood$Density[g]
            fluxes_wp$Moisture[t]<- wood$Moisture[g]
            fluxes_wp$Length_S[t]<- wood$Length_S[g]
            fluxes_wp$Moisture_H[t]<- wood$Moisture_H[g]
            fluxes_wp$Moisture_S[t]<- wood$Moisture_S[g]
          }
        }
      }
      
      lmm_co2_lme_150<- lme(log(CO2_flux)~Temp+VWC+Diam+Moisture+Density+Length_S+Moisture_S,
                            random=~1|Tree,data=fluxes_wp)
  #CO2 model with wood properties
      # model with temporal autocorrelation
      lmm_co2_AR1_150 <- update(lmm_co2_lme_150,cor=corAR1(form=~DOY|Tree))

      acf(residuals(lmm_co2_AR1_150,type='normalized'))
      
      plot(lmm_co2_AR1_150)
      summary(lmm_co2_AR1_150)
      MuMIn::r.squaredGLMM(lmm_co2_AR1_150)
      car::Anova(lmm_co2_AR1_150,3)
      
  #CH4 model with wood properties    
      # use box-cox transformed ch4 fluxt
      lmm_ch4bc_lme_150<- lme(CH4_flux_bc~Temp+VWC+Diam+Moisture+Density+Length_S+Moisture_S,
                          random=~1|Tree,data=fluxes_wp)
      
      # add corAR1
      lmm_ch4bc_AR1_150 <- update(lmm_ch4bc_lme_150,
                              cor=corAR1(form=~DOY|Tree))
      plot(lmm_ch4bc_AR1_150)
      anova(lmm_ch4bc_lme_150,lmm_ch4bc_AR1_150)

      acf(residuals(lmm_ch4bc_AR1_150,type='normalized'))
      
      plot(lmm_ch4bc_AR1_150)
      summary(lmm_ch4bc_AR1_150)
      MuMIn::r.squaredGLMM(lmm_ch4bc_AR1_150)
      car::Anova(lmm_ch4bc_AR1_150,3)
      
# Correlation between manual CO2 and CH4 fluxes
      fluxes<-read.table("Manual_fluxes.txt", fill=TRUE,header=TRUE,sep="\t",na.strings=c("NA","#N/A!","#N/A","#NA!"))
      
      CH4_flux_tree_fixed<-ggplot(fluxes, aes(x=log(CO2_flux), y=CH4_flux)) +
        #geom_line() +
        geom_point() +
        geom_smooth(method='glm', span=0.1, size=1.5, na.rm=T) +
        #scale_color_brewer(palette="Dark2", labels = c("50cm","100cm","150cm"), name="Stem height") +  #Colors displayed for this palette: #1b9e77 for green, #d95f02 for orange and #7570b3 for purple
        #scale_shape(labels=c("50cm","100cm","150cm"), name="Stem height") +
        #geom_hline(yintercept = 0, linetype=2, size=0.5) +
        facet_wrap(~Tree, scales="free") +
        xlab("DOY") + 
        ylab(bquote('CH'[4]~ 'flux  (nmol ' ~CH[4]~ m^-2~s^-1*')')) +
        theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
              axis.text.y = element_text(colour = "grey20", size = 12),
              text = element_text(size = 16)) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        theme(legend.position = c(1, 0),
              legend.justification = c(1.5, -0.2),
              legend.title = element_text(face =2),
              legend.text=element_text(size=12),
              axis.title=element_text(size=12),
              axis.text = element_text(color = "black", size=10)) +
        guides(color = guide_legend(reverse = TRUE), shape = guide_legend(reverse = TRUE)) 

      
      ggplot(df, aes(x= new_price, y= carat, color = cut)) +
        geom_point(alpha = 0.3) +
        facet_wrap(~clarity, scales = "free_y") +
        geom_smooth(method = "lm", formula = formula, se = F) +
        stat_poly_eq(aes(label = paste(..rr.label..)), 
                     label.x.npc = "right", label.y.npc = 0.15,
                     formula = formula, parse = TRUE, size = 3)+
        stat_fit_glance(method = 'lm',
                        method.args = list(formula = formula),
                        geom = 'text',
                        aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),
                        label.x.npc = 'right', label.y.npc = 0.35, size = 3)
      
    
      CH4_flux_tree_fixed<-ggplot(fluxes, aes(x=(CO2_flux), y=CH4_flux)) +
        geom_point() +
        #xlim(0,7) +
        ylim(-1,5)
      
      Tree_A<-lm(CH4_flux~CH4_flux_bc_z, data=fluxes)
      summary(Tree_A)$adj.r.squared
      summary(Tree_A)
      
      Tree_A<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="A"))
      summary(Tree_A)$adj.r.squared
      summary(Tree_A)
      
      Tree_B<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="B"))
      summary(Tree_B)$adj.r.squared
      summary(Tree_B)

      Tree_C<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="C"))
      summary(Tree_C)$adj.r.squared
      summary(Tree_C)
      
      Tree_D<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="D"))
      summary(Tree_D)$adj.r.squared
      summary(Tree_D)
      
      Tree_E<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="E"))
      summary(Tree_E)$adj.r.squared
      summary(Tree_E)
      
      Tree_F<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="F"))
      summary(Tree_F)$adj.r.squared
      summary(Tree_F)
      
      Tree_G<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="G"))
      summary(Tree_G)$adj.r.squared
      summary(Tree_G)
      
      Tree_H<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="H"))
      summary(Tree_H)$adj.r.squared
      summary(Tree_H)
      
      Tree_I<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="I"))
      summary(Tree_I)$adj.r.squared
      summary(Tree_I)
      
      Tree_J<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="J"))
      summary(Tree_J)$adj.r.squared
      summary(Tree_J)
      
      Tree_K<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="K"))
      summary(Tree_K)$adj.r.squared
      summary(Tree_K)
      
      Tree_L<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="L"))
      summary(Tree_L)$adj.r.squared
      summary(Tree_L)
      
      Tree_M<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="M"))
      summary(Tree_M)$adj.r.squared
      summary(Tree_M)
      
      Tree_N<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="N"))
      summary(Tree_N)$adj.r.squared
      summary(Tree_N)
      
      Tree_O<-lm(CH4_flux~CO2_flux, data=subset(fluxes,Tree=="O"))
      summary(Tree_O)$adj.r.squared
      summary(Tree_O)
 
# EXPLORATORY ANALYSES AND PLOTS ----   
    
    # CO2-CH4 relation using all the manual measurements  
    CO2_CH4<-ggplot(fluxes, aes(x=CO2_flux, y=(log(7+CH4_flux)))) + #lof(7+CH4_flux) to avoid negative values. This should be adjusted once I have the final maual fluxes
      geom_point(na.rm=TRUE, size=1, alpha=0.7) +
      #facet_wrap(~Tree, scales="free") +
      #scale_color_manual(values=c("red4", "dodgerblue4", "black")) + 
      stat_smooth(method = 'lm', na.rm=TRUE, se=F) +
      xlab('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') + 
      ylab('CH'[4]~ 'flux  (nmol '~CH[4]~ m^-2~s^-1*')') +
      #annotate("text", x=1.8, y=0.9, label= "a)",size=4.5, fontface =2) +
      #annotate("text", x=13, y=0.8, label= "UpperStem",color="black",size=4, fontface =2) +
      #annotate("text", x=13, y=0.65, label= "LowerStem",color="dodgerblue4",size=4, fontface =2) +
      #annotate("text", x=13, y=0.5, label= "Soil",color="red4",size=4, fontface =2) +
      #annotate("text", x=13, y=-0.4, label= CO2.CH4.soil,color="red4",size=4, parse=TRUE) +
      #annotate("text", x=8.7, y=0.7, label= CO2.CH4.lower,color="dodgerblue4",size=4, parse=TRUE) +
      #annotate("text", x=4, y=0, label= CO2.CH4.upper,color="black",size=4, parse=TRUE) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=11),
            axis.text = element_text(color = "black"))
    
    sumfun <- function(x, ...){
      c(m=mean(x, ...), sd=sd(x, ...), l=length(x))
    }
    daily_manual <- summaryBy(CO2_flux + CH4_flux ~ Tree + Height + Diam + Location, data=fluxes,FUN=c(sumfun), na.rm=TRUE)
    daily_manual$CH4.CI<-qnorm(0.975)*daily_manual$CH4_flux.sd/sqrt(daily_manual$CH4_flux.l)
    daily_manual$CO2.CI<-qnorm(0.975)*daily_manual$CO2_flux.sd/sqrt(daily_manual$CO2_flux.l)
    
    # CO2-CH4 relation using annual means for each tree and height. Similar results when using just annual means by tree  
    CO2_CH4<-ggplot(daily_manual, aes(x=CO2_flux.m, y=(log(1+CH4_flux.m)))) + #lof(7+CH4_flux) to avoid negative values. This should be adjusted once I have the final maual fluxes
      geom_point(na.rm=TRUE, size=1, alpha=0.7) +
      #facet_wrap(~Tree, scales="free") +
      #scale_color_manual(values=c("red4", "dodgerblue4", "black")) + 
      stat_smooth(method = 'lm', na.rm=TRUE, se=F) +
      xlab('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') + 
      ylab('CH'[4]~ 'flux  (nmol '~CH[4]~ m^-2~s^-1*')') +
      #annotate("text", x=1.8, y=0.9, label= "a)",size=4.5, fontface =2) +
      #annotate("text", x=13, y=0.8, label= "UpperStem",color="black",size=4, fontface =2) +
      #annotate("text", x=13, y=0.65, label= "LowerStem",color="dodgerblue4",size=4, fontface =2) +
      #annotate("text", x=13, y=0.5, label= "Soil",color="red4",size=4, fontface =2) +
      #annotate("text", x=13, y=-0.4, label= CO2.CH4.soil,color="red4",size=4, parse=TRUE) +
      #annotate("text", x=8.7, y=0.7, label= CO2.CH4.lower,color="dodgerblue4",size=4, parse=TRUE) +
      #annotate("text", x=4, y=0, label= CO2.CH4.upper,color="black",size=4, parse=TRUE) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=11),
            axis.text = element_text(color = "black"))
    
    # CO2 and CH4 by stem height ----
    #Relation between CH4 emissions and stem heigh using all the manual data. Facet by tree
    CH4<-ggplot(fluxes, aes(x=Height, y=CH4_flux)) +
      geom_point(na.rm=TRUE, size=1, alpha=0.7) +
      stat_smooth(method = 'nls', formula = y ~ exp(a *x+b),
                  method.args = list(start=c(a=0.1,b=1)), na.rm=TRUE, se=F) +
      coord_flip() +
      #scale_x_continuous(limits = c(20,180),breaks = c(50, 100, 150), label = c("50", "100", "150")) +
      
      facet_wrap(~Tree, scales="free") +
      #scale_color_manual(values=c("red4", "dodgerblue4", "black")) + 
      #stat_smooth(method = 'lm', na.rm=TRUE, se=F) +
      #xlab('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') + 
      ylab('CH'[4]~ 'flux  (nmol '~CH[4]~ m^-2~s^-1*')') +
      #annotate("text", x=1.8, y=0.9, label= "a)",size=4.5, fontface =2) +
      #annotate("text", x=13, y=0.8, label= "UpperStem",color="black",size=4, fontface =2) +
      #annotate("text", x=13, y=0.65, label= "LowerStem",color="dodgerblue4",size=4, fontface =2) +
      #annotate("text", x=13, y=0.5, label= "Soil",color="red4",size=4, fontface =2) +
      #annotate("text", x=13, y=-0.4, label= CO2.CH4.soil,color="red4",size=4, parse=TRUE) +
      #annotate("text", x=8.7, y=0.7, label= CO2.CH4.lower,color="dodgerblue4",size=4, parse=TRUE) +
      #annotate("text", x=4, y=0, label= CO2.CH4.upper,color="black",size=4, parse=TRUE) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=11),
            axis.text = element_text(color = "black"))
    
    #Relation between CH4 emissions and stem height using the mean for each tree and plotting the confidence interval from all measurements within tree
    
    CH4<-ggplot(subset(daily_manual,CH4_flux.m<10), aes(x=Height, y=CH4_flux.m)) +
      geom_point(na.rm=TRUE, size=1, alpha=0.7) +
      geom_linerange(aes(ymin = CH4_flux.m-CH4.CI, ymax = CH4_flux.m+CH4.CI)) +
      #Eliminar les dues linies seguents si la relacio no es significativa    
      stat_smooth(method = 'nls', formula = y ~ exp(a *x+b),
                  method.args = list(start=c(a=0.1,b=1)), na.rm=TRUE, se=F) +
      
      #geom_errorbar(x=Height,ymin=daily_manual$CH4_flux.m-daily_manual$CH4.CI,ymax=daily_manual$CH4_flux.m+daily_manual$CH4.CI,position="dodge",width = 1.8),size=0.5) +     
      #stat_smooth(method = 'nls', formula = y ~ exp(a *x+b),
      #            method.args = list(start=c(a=0.1,b=1)), na.rm=TRUE, se=F) +
      coord_flip() +
      #scale_x_continuous(limits = c(20,180),breaks = c(50, 100, 150), label = c("50", "100", "150")) +
      
      #facet_wrap(~Tree, scales="free") +
      #scale_color_manual(values=c("red4", "dodgerblue4", "black")) + 
      stat_smooth(method = 'lm', na.rm=TRUE, se=F) +
      #xlab('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') + 
      ylab('CH'[4]~ 'flux  (nmol '~CH[4]~ m^-2~s^-1*')') +
      #annotate("text", x=1.8, y=0.9, label= "a)",size=4.5, fontface =2) +
      #annotate("text", x=13, y=0.8, label= "UpperStem",color="black",size=4, fontface =2) +
      #annotate("text", x=13, y=0.65, label= "LowerStem",color="dodgerblue4",size=4, fontface =2) +
      #annotate("text", x=13, y=0.5, label= "Soil",color="red4",size=4, fontface =2) +
      #annotate("text", x=13, y=-0.4, label= CO2.CH4.soil,color="red4",size=4, parse=TRUE) +
      #annotate("text", x=8.7, y=0.7, label= CO2.CH4.lower,color="dodgerblue4",size=4, parse=TRUE) +
      #annotate("text", x=4, y=0, label= CO2.CH4.upper,color="black",size=4, parse=TRUE) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=11),
            axis.text = element_text(color = "black"))
    
    
    hola<-lm(CH4_flux.m~Height,daily_manual)
    hola<-lm(CH4_flux.m~Height,subset(daily_manual,CH4_flux.m<10))
    
    #Relation between CO2 emissions and stem heigh using all the manual data
    CO2<-ggplot(fluxes, aes(x=Height, y=CO2_flux)) +
      geom_point(na.rm=TRUE, size=1, alpha=0.7) +
      stat_smooth(method = 'nls', formula = y ~ exp(a *x+b),
                  method.args = list(start=c(a=0.1,b=1)), na.rm=TRUE, se=F) +
      coord_flip() +
      #scale_x_continuous(limits = c(20,180),breaks = c(50, 100, 150), label = c("50", "100", "150")) +
      
      #facet_wrap(~Tree, scales="free") +
      #scale_color_manual(values=c("red4", "dodgerblue4", "black")) + 
      #stat_smooth(method = 'lm', na.rm=TRUE, se=F) +
      ylab('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') + 
      #ylab('CH'[4]~ 'flux  (nmol '~CH[4]~ m^-2~s^-1*')') +
      #annotate("text", x=1.8, y=0.9, label= "a)",size=4.5, fontface =2) +
      #annotate("text", x=13, y=0.8, label= "UpperStem",color="black",size=4, fontface =2) +
      #annotate("text", x=13, y=0.65, label= "LowerStem",color="dodgerblue4",size=4, fontface =2) +
      #annotate("text", x=13, y=0.5, label= "Soil",color="red4",size=4, fontface =2) +
      #annotate("text", x=13, y=-0.4, label= CO2.CH4.soil,color="red4",size=4, parse=TRUE) +
      #annotate("text", x=8.7, y=0.7, label= CO2.CH4.lower,color="dodgerblue4",size=4, parse=TRUE) +
      #annotate("text", x=4, y=0, label= CO2.CH4.upper,color="black",size=4, parse=TRUE) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=11),
            axis.text = element_text(color = "black"))
    
    #Relation between CO2 emissions and stem height using the mean for each tree and plotting the confidence interval from all measurements within tree
    
    CO2<-ggplot(daily_manual, aes(x=Height, y=CO2_flux.m)) +
      geom_point(na.rm=TRUE, size=1, alpha=0.7) +
      geom_linerange(aes(ymin = CO2_flux.m-CO2.CI, ymax = CO2_flux.m+CO2.CI)) +
      #Eliminar les dues linies seguents si la relacio no es significativa    
      stat_smooth(method = 'nls', formula = y ~ exp(a *x+b),
                  method.args = list(start=c(a=0.1,b=1)), na.rm=TRUE, se=F) +
      
      #geom_errorbar(x=Height,ymin=daily_manual$CH4_flux.m-daily_manual$CH4.CI,ymax=daily_manual$CH4_flux.m+daily_manual$CH4.CI,position="dodge",width = 1.8),size=0.5) +     
      #stat_smooth(method = 'nls', formula = y ~ exp(a *x+b),
      #            method.args = list(start=c(a=0.1,b=1)), na.rm=TRUE, se=F) +
      coord_flip() +
      #scale_x_continuous(limits = c(20,180),breaks = c(50, 100, 150), label = c("50", "100", "150")) +
      
      facet_wrap(~Tree, scales="free") +
      #scale_color_manual(values=c("red4", "dodgerblue4", "black")) + 
      #stat_smooth(method = 'lm', na.rm=TRUE, se=F) +
      ylab('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') + 
      #ylab('CH'[4]~ 'flux  (nmol '~CH[4]~ m^-2~s^-1*')') +
      #annotate("text", x=1.8, y=0.9, label= "a)",size=4.5, fontface =2) +
      #annotate("text", x=13, y=0.8, label= "UpperStem",color="black",size=4, fontface =2) +
      #annotate("text", x=13, y=0.65, label= "LowerStem",color="dodgerblue4",size=4, fontface =2) +
      #annotate("text", x=13, y=0.5, label= "Soil",color="red4",size=4, fontface =2) +
      #annotate("text", x=13, y=-0.4, label= CO2.CH4.soil,color="red4",size=4, parse=TRUE) +
      #annotate("text", x=8.7, y=0.7, label= CO2.CH4.lower,color="dodgerblue4",size=4, parse=TRUE) +
      #annotate("text", x=4, y=0, label= CO2.CH4.upper,color="black",size=4, parse=TRUE) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=11),
            axis.text = element_text(color = "black"))
    
    
    hola<-lm(CH4_flux.m~Height,daily_manual)
    hola<-lm(CH4_flux.m~Height,subset(daily_manual,CH4_flux.m<10))
    
    
    # CH4 and CO2 with diameter
    #Relation between CH4 emissions and stem heigh using all the manual data
    CH4<-ggplot(subset(fluxes,CH4_flux<80), aes(x=Diam, y=CH4_flux)) +
      geom_point(na.rm=TRUE, size=1, alpha=0.7) +
      stat_smooth(method = 'nls', formula = y ~ exp(a *x+b),
                  method.args = list(start=c(a=0.1,b=1)), na.rm=TRUE, se=F) +
      #coord_flip() +
      #scale_x_continuous(limits = c(20,180),breaks = c(50, 100, 150), label = c("50", "100", "150")) +
      
      #facet_wrap(~Location, scales="free") +
      #scale_color_manual(values=c("red4", "dodgerblue4", "black")) + 
      #stat_smooth(method = 'lm', na.rm=TRUE, se=F) +
      #xlab('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') + 
      ylab('CH'[4]~ 'flux  (nmol '~CH[4]~ m^-2~s^-1*')') +
      #annotate("text", x=1.8, y=0.9, label= "a)",size=4.5, fontface =2) +
      #annotate("text", x=13, y=0.8, label= "UpperStem",color="black",size=4, fontface =2) +
      #annotate("text", x=13, y=0.65, label= "LowerStem",color="dodgerblue4",size=4, fontface =2) +
      #annotate("text", x=13, y=0.5, label= "Soil",color="red4",size=4, fontface =2) +
      #annotate("text", x=13, y=-0.4, label= CO2.CH4.soil,color="red4",size=4, parse=TRUE) +
      #annotate("text", x=8.7, y=0.7, label= CO2.CH4.lower,color="dodgerblue4",size=4, parse=TRUE) +
      #annotate("text", x=4, y=0, label= CO2.CH4.upper,color="black",size=4, parse=TRUE) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=11),
            axis.text = element_text(color = "black"))
    
    #Relation between CH4 emissions and stem height using the mean for each tree and plotting the confidence interval from all measurements within tree
    
    CH4<-ggplot(subset(daily_manual,CH4_flux.m<10), aes(x=Diam, y=CH4_flux.m)) +
      geom_point(na.rm=TRUE, size=1, alpha=0.7) +
      geom_linerange(aes(ymin = CH4_flux.m-CH4.CI, ymax = CH4_flux.m+CH4.CI)) +
      #Eliminar les dues linies seguents si la relacio no es significativa    
      stat_smooth(method = 'nls', formula = y ~ exp(a *x+b),
                  method.args = list(start=c(a=0.1,b=1)), na.rm=TRUE, se=F) +
      
      #geom_errorbar(x=Height,ymin=daily_manual$CH4_flux.m-daily_manual$CH4.CI,ymax=daily_manual$CH4_flux.m+daily_manual$CH4.CI,position="dodge",width = 1.8),size=0.5) +     
      #stat_smooth(method = 'nls', formula = y ~ exp(a *x+b),
      #            method.args = list(start=c(a=0.1,b=1)), na.rm=TRUE, se=F) +
      #coord_flip() +
      #scale_x_continuous(limits = c(20,180),breaks = c(50, 100, 150), label = c("50", "100", "150")) +
      
      #facet_wrap(~Tree, scales="free") +
      #scale_color_manual(values=c("red4", "dodgerblue4", "black")) + 
      stat_smooth(method = 'lm', na.rm=TRUE, se=F) +
      #xlab('CO'[2]~ 'flux  ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')') + 
      #ylab('CH'[4]~ 'flux  (nmol '~CH[4]~ m^-2~s^-1*')') +
      #annotate("text", x=1.8, y=0.9, label= "a)",size=4.5, fontface =2) +
      #annotate("text", x=13, y=0.8, label= "UpperStem",color="black",size=4, fontface =2) +
      #annotate("text", x=13, y=0.65, label= "LowerStem",color="dodgerblue4",size=4, fontface =2) +
      #annotate("text", x=13, y=0.5, label= "Soil",color="red4",size=4, fontface =2) +
      #annotate("text", x=13, y=-0.4, label= CO2.CH4.soil,color="red4",size=4, parse=TRUE) +
      #annotate("text", x=8.7, y=0.7, label= CO2.CH4.lower,color="dodgerblue4",size=4, parse=TRUE) +
      #annotate("text", x=4, y=0, label= CO2.CH4.upper,color="black",size=4, parse=TRUE) +
      theme_bw() + theme(panel.grid = element_blank()) +
      theme(legend.position="none",
            axis.title=element_text(size=11),
            axis.text = element_text(color = "black"))

    
