library(tidyverse)
library(ggplot2)
library(latex2exp)
library(broom) #used for the tidy function
library(xtable)
library(gridExtra)
library(grid)
library(gtable)
library(magick)
setwd("~/Desktop/sandiaProjects/languageModels/simulations")


ytitle = TeX("Monte Carlo Estimate of Mean $\\| \\cdot \\|_1$ Error")
find_cell <- function(table,row,col,name="core-bg") {
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}

sim2 <- function() {
  dat = read_csv("data/sim2.csv")
  numMC = dat$M[1]
  #quick plot to assess noise level
  #ggplot(data = dat,aes(x = mcVal))+geom_histogram()+facet_wrap(~n)
  
  xlabel = TeX("$\\log(n)$")
  #xlabel = ""
  #ylabel = ytitle
  ylabel = ""
  
  zVal = qnorm(.975)
  averaged = dat %>%
    group_by(beta,s,n,type) %>%
    summarise(expL1 = mean(mcVal), sdL1 = sd(mcVal),tolerance = (zVal/sqrt(numMC))*(sdL1/expL1)) %>%
    mutate(logn = log(n),logy = log(expL1),lower = log(expL1)-tolerance,upper = log(expL1)+tolerance,weight = (expL1^2)/(sdL1^2)) %>%
    mutate(beta_s = str_c(beta,s,sep = ",")) %>%
    mutate(s = as.factor(s)) %>%
    mutate(Estimator = type) %>%
    filter(s != 4)
  mcPlot = ggplot(data = averaged,aes(x = logn,y = logy,ymin = lower,ymax = upper,color=Estimator))+geom_point(size=.7)+geom_line(data=averaged,aes(linetype = s))+geom_errorbar()+
    xlab(xlabel)+ylab(ylabel)+
    ggtitle("Simulation 2: High s and \U03B2=1")+
    scale_x_continuous(limits=c(4.5,14.5),breaks = 5:14)+
    scale_y_continuous(limits=c(-12,-2),breaks = seq(-12,-2,by=2))+
    scale_color_manual(guide = "none",labels = c("EPE","SS","TSS"),values = c("#7CAE00","#F8766D","#00BFC4"))+
    theme(
      #LABLES APPEARANCE
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      plot.title = element_text(hjust = 0.5, size=43, face= "bold", colour= "black" ),
      axis.title.x = element_text(size=28, face="bold", colour = "black"),    
      axis.title.y = element_text(size=28, face="bold", colour = "black"),    
      axis.text.x = element_text(size=21, face="bold", colour = "black"), 
      axis.text.y = element_text(size=21, face="bold", colour = "black"),
      strip.text.x = element_text(size = 11, face="bold", colour = "black" ),
      strip.text.y = element_text(size = 11, face="bold", colour = "black"),
      axis.line.x = element_line(color="black", size = 0.3),
      axis.line.y = element_line(color="black", size = 0.3),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      legend.title = element_text(size = 43),
      legend.text = element_text(size = 33),
      legend.position = "bottom"
    )+
    guides(linetype = guide_legend(override.aes = list(linewidth=1.2)))
    
  
  #Now we create the table
  slopeTable = averaged %>%
    group_by(type,s) %>%
    do(tidy(lm(logy ~ logn,data=.))) %>%
    mutate(estimate = round(estimate,3)) %>%
    mutate(lower = estimate - qnorm(.975)*std.error, upper = estimate + qnorm(.975)*std.error) %>%
    filter(term == "logn") %>%
    rename(Estimator = type,'OLS_Slope_Estimate' = estimate) %>%
    select(Estimator,s,OLS_Slope_Estimate) %>%
    arrange(s) %>% 
    tableGrob(rows=NULL) #%>%
    #select(Estimator,s,estimate,lower,upper) %>%
    #rename("lower 95% Confidence Bound" = lower,"upper 95% Confidence Bound" = upper) %>%
    #select(Estimator,s,estimate) %>%
    #xtable()
  
  #change colors to easier indicate what slopes match to which lines
  ind1 = find_cell(slopeTable,2,1)
  ind2 = find_cell(slopeTable,5,1)
  slopeTable$grobs[ind1][[1]][["gp"]] <- gpar(fill = "#7CAE00",lwd=5)
  slopeTable$grobs[ind2][[1]][["gp"]] <- gpar(fill = "#7CAE00",lwd=5)
  ind3 = find_cell(slopeTable,3,1)
  ind4 = find_cell(slopeTable,6,1)
  slopeTable$grobs[ind3][[1]][["gp"]] <- gpar(fill = "#F8766D",lwd=5)
  slopeTable$grobs[ind4][[1]][["gp"]] <- gpar(fill = "#F8766D",lwd=5)
  ind5 = find_cell(slopeTable,4,1)
  ind6 = find_cell(slopeTable,7,1)
  slopeTable$grobs[ind5][[1]][["gp"]] <- gpar(fill = "#00BFC4",lwd=5)
  slopeTable$grobs[ind6][[1]][["gp"]] <- gpar(fill = "#00BFC4",lwd=5)
  
  #mcPlot = mcPlot+annotation_custom(slopeTable,xmin=11.8,xmax=15,ymin=-3.2,ymax=-2.2)
  mcPlot = mcPlot+annotation_custom(slopeTable,xmin=12.1,xmax=15,ymin=-3.1,ymax=-2.4)
  return (mcPlot)

}


sim3 <- function() {
  dat = read_csv("data/sim3.csv")
  numMC = dat$M[1]
  #quick plot to assess noise level
  #ggplot(data = dat,aes(x = mcVal))+geom_histogram()+facet_wrap(~n)
  
  xlabel = TeX("$\\log(n)$")
  #ylabel = ytitle
  ylabel = ""
  
  zVal = qnorm(.975)
  averaged = dat %>%
    group_by(beta,s,n,type) %>%
    summarise(expL1 = mean(mcVal), sdL1 = sd(mcVal),tolerance = (zVal/sqrt(numMC))*(sdL1/expL1)) %>%
    mutate(logn = log(n),logy = log(expL1),lower = log(expL1)-tolerance,upper = log(expL1)+tolerance,weight = (expL1^2)/(sdL1^2)) %>%
    mutate(beta_s = str_c(beta,s,sep = ",")) %>%
  
      mutate(s = as.factor(s)) %>%
    mutate(Estimator = type)
  mcPlot = ggplot(data = averaged,aes(x = logn,y = logy,ymin = lower,ymax = upper,color=Estimator))+geom_point(size=.7)+geom_line(data=averaged,aes(linetype=s))+geom_errorbar()+
    xlab(xlabel)+ylab(ylabel)+
    ggtitle("Simulation 3: s=1.5 and \U03B2=1/s")+
    scale_x_continuous(limits=c(4.5,14.5),breaks = 5:14)+
    scale_y_continuous(limits=c(-4.5,-1.1),breaks = c(-4,-3,-2))+
    scale_color_manual(guide = "none",labels = c("EPE","SS","TSS"),values = c("#7CAE00","#F8766D","#00BFC4"))+
    theme(
      #LABLES APPEARANCE
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      plot.title = element_text(hjust = 0.5, size=43, face= "bold", colour= "black" ),
      axis.title.x = element_text(size=28, face="bold", colour = "black"),    
      axis.title.y = element_text(size=28, face="bold", colour = "black"),    
      axis.text.x = element_text(size=21, face="bold", colour = "black"), 
      axis.text.y = element_text(size=21, face="bold", colour = "black"),
      strip.text.x = element_text(size = 11, face="bold", colour = "black" ),
      strip.text.y = element_text(size = 11, face="bold", colour = "black"),
      axis.line.x = element_line(color="black", size = 0.3),
      axis.line.y = element_line(color="black", size = 0.3),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      legend.title = element_text(size = 43),
      legend.text = element_text(size = 33),
      legend.position = "bottom"
    )+
    guides(linetype = guide_legend(override.aes = list(linewidth=1.2)))
    
  
  #Now we create the table
  slopeTable = averaged %>%
    group_by(type) %>%
    do(tidy(lm(logy ~ logn,data=.))) %>%
    mutate(estimate = round(estimate,digits = 3)) %>%
    mutate(lower = estimate - qnorm(.975)*std.error, upper = estimate + qnorm(.975)*std.error) %>%
    filter(term == "logn") %>%
    rename(Estimator = type,'OLS_Slope_Estimate' = estimate) %>%
    select(Estimator,OLS_Slope_Estimate) %>% 
    tableGrob(rows=NULL) #%>%
    #select(Estimator,s,estimate,lower,upper) %>%
    #rename("lower 95% Confidence Bound" = lower,"upper 95% Confidence Bound" = upper) %>%
    #select(Estimator,s,estimate) %>%
    #xtable()
  
  #change colors to easier indicate what slopes match to which lines
  ind1 = find_cell(slopeTable,2,1)
  slopeTable$grobs[ind1][[1]][["gp"]] <- gpar(fill = "#7CAE00",lwd=5)
  ind2 = find_cell(slopeTable,3,1)
  slopeTable$grobs[ind2][[1]][["gp"]] <- gpar(fill = "#F8766D",lwd=5)
  ind3 = find_cell(slopeTable,4,1)
  slopeTable$grobs[ind3][[1]][["gp"]] <- gpar(fill = "#00BFC4",lwd=5)
  
  
  mcPlot = mcPlot+annotation_custom(slopeTable,xmin=12.88,xmax=14.5,ymin=-1.45,ymax=-.90)
  return (mcPlot)
}
#width = 900
#height = 800


sim1 <- function() {
  dat = read_csv("data/sim1.csv")
  numMC = dat$M[1]
  #quick plot to assess noise level
  #ggplot(data = dat,aes(x = mcVal))+geom_histogram()+facet_wrap(~n)
  
  xlabel = TeX("$\\log(n)$")
  #xlabel = ""
  ylabel = ytitle
  
  zVal = qnorm(.975)
  averaged = dat %>%
    group_by(beta,s,n,type) %>%
    summarise(expL1 = mean(mcVal), sdL1 = sd(mcVal),tolerance = (zVal/sqrt(numMC))*(sdL1/expL1)) %>%
    mutate(logn = log(n),logy = log(expL1),lower = log(expL1)-tolerance,upper = log(expL1)+tolerance,weight = (expL1^2)/(sdL1^2)) %>%
    mutate(beta_s = str_c(beta,s,sep = ",")) %>%
    mutate(s = as.factor(s)) %>%
    mutate(Estimator = type)
  mcPlot = ggplot(data = averaged,aes(x = logn,y = logy,ymin = lower,ymax = upper,color=Estimator))+geom_point(size=.7)+geom_line(data=averaged,aes(linetype=s))+geom_errorbar()+
    xlab(xlabel)+ylab(ylabel)+
    ggtitle("Simulation 1: s=1.05 and \U03B2=1/(s+3))")+
    scale_x_continuous(limits=c(4.5,14.5),breaks = 5:14)+
    scale_y_continuous(limits=c(-12.5,-2.3),breaks=seq(-12,-2,by=2))+
    scale_color_manual(guide="none",labels = c("EPE","SS"),values = c("#7CAE00","#F8766D"))+
    theme(
      #LABLES APPEARANCE
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      plot.title = element_text(hjust = 0.5, size=43, face= "bold", colour= "black" ),
      axis.title.x = element_text(size=28, face="bold", colour = "black"),    
      axis.title.y = element_text(size=28, face="bold", colour = "black"),    
      axis.text.x = element_text(size=21, face="bold", colour = "black"), 
      axis.text.y = element_text(size=21, face="bold", colour = "black"),
      strip.text.x = element_text(size = 11, face="bold", colour = "black" ),
      strip.text.y = element_text(size = 11, face="bold", colour = "black"),
      axis.line.x = element_line(color="black", size = 0.3),
      axis.line.y = element_line(color="black", size = 0.3),
      panel.border = element_rect(colour = "black", fill=NA, size=0.3),
      legend.title = element_text(size = 43),
      legend.text = element_text(size = 33),
      legend.position = "bottom"
    )+
    guides(linetype = guide_legend(override.aes = list(linewidth=1.2)))
    
  
  ColorTable <- data.frame(Estimator = c("EPE","SS","TSS")) %>% tableGrob(rows=NULL)
  #Now we create the table
  #slopeTable = averaged %>%
  #  group_by(type) %>%
  #  do(tidy(lm(logy ~ logn,data=.))) %>%
  #  mutate(estimate = round(estimate,digits = 3)) %>%
  #  mutate(lower = estimate - qnorm(.975)*std.error, upper = estimate + qnorm(.975)*std.error) %>%
  #  filter(term == "logn") %>%
  #  rename(Estimator = type,'OLS_Slope_Estimate' = estimate) %>%
  #  select(Estimator,OLS_Slope_Estimate) %>% 
  #  tableGrob(rows=NULL) #%>%
  #  #select(Estimator,s,estimate,lower,upper) %>%
  #  #rename("lower 95% Confidence Bound" = lower,"upper 95% Confidence Bound" = upper) %>%
  #  #select(Estimator,s,estimate) %>%
  #  #xtable()
  #
  #change colors to easier indicate what slopes match to which lines
  ind1 = find_cell(ColorTable,2,1)
  ColorTable$grobs[ind1][[1]][["gp"]] <- gpar(fill = "#7CAE00",lwd=5)
  ind2 = find_cell(ColorTable,3,1)
  ColorTable$grobs[ind2][[1]][["gp"]] <- gpar(fill = "#F8766D",lwd=5)
  ind3 = find_cell(ColorTable,4,1)
  ColorTable$grobs[ind3][[1]][["gp"]] <- gpar(fill = "#00BFC4",lwd=5)
  
  
  mcPlot = mcPlot+annotation_custom(ColorTable,xmin=14.075,xmax=15,ymin=-2.7,ymax=-2.34)
  return (mcPlot)
}


#WLS slopes -- the weights are calculated as the inverses of the CLT+delta method variances.
#The CLT delta method variances are sigma_n^2 / mu_n^2. These are estimated using data.
#I don't really know how reliable those standard errors being produced are by the lm function.


png("data/sim2.png",width = 900,height = 800,units = "px")
p2 = sim2()
p2
dev.off()
###########
png("data/sim3.png",width = 900, height = 800,units = "px")
p3 = sim3()
p3
dev.off()
##########
png("data/sim1.png",width = 900, height = 800,units = "px")
p1 = sim1()
p1
dev.off()

#################
#Image Concatenation
p1 = image_read("data/sim1.png")
p2 = image_read("data/sim2.png")
p3 = image_read("data/sim3.png")
newImg = image_append(c(p1,p2,p3),stack = FALSE)
image_write(newImg,path = "data/mcSims.png",format="png")






