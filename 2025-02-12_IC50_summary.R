library(tidyverse)
library(readxl)
library(ggsci)
library(growthcurver)
library(MESS)
library(drc)


#metadata
##########
read_metadata <- function(metadata_file){
  #jmena jamek
  well <- paste(rep(LETTERS[1:8], 12), 
                rep(c(1:12), each = 8),
                sep = "")
  md <- read_excel(metadata_file, col_names = FALSE, range = cell_cols("A:M"))
  
  ###precursor
  start_row <- which(md == "<>precursor")
  excel_range <- paste("B", start_row, ":M", start_row + 8, sep = "")
  precursor <- read_excel(metadata_file, range = excel_range, na = "NA")
  precursor <- as.vector(as.matrix(precursor))
  
  ###derivate
  start_row <- which(md == "<>derivate")
  excel_range <- paste("B", start_row, ":M", start_row + 8, sep = "")
  derivate <- read_excel(metadata_file, range = excel_range, na = "NA")
  derivate <- as.vector(as.matrix(derivate))
  
  ###concentrations
  start_row <- which(md == "<>concentration")
  excel_range <- paste("B", start_row, ":M", start_row + 8, sep = "")
  concentrations <- read_excel(metadata_file, range = excel_range, na = "NA")
  concentrations <- as.vector(as.matrix(concentrations))
  concentrations <- as.numeric(concentrations)
  
  ###strains
  start_row <- which(md == "<>strain")
  excel_range <- paste("B", start_row, ":M", start_row + 8, sep = "")
  strains <- read_excel(metadata_file, range = excel_range, na = "NA")
  strains <- as.vector(as.matrix(strains))
  
  #data_frame
  metadata <- data.frame(well, precursor, derivate, concentrations, strains)
  names(metadata) <- c("well", "precursor", "derivate", "concentration", "strain")
  
  return(metadata)
}

metadata.nema <- read_metadata('2025-01-14_Xnema/2025-01-14_plate_metadata.xlsx')
metadata.bovi <- read_metadata('2025-02-06_Xbovi02/2025-02-06_plate_metadata.xlsx')
metadata.plau <- read_metadata('2025-01-31_Plaumo/2025-01-31_plate_metadata.xlsx')

metadata <- rbind(metadata.bovi, metadata.nema, metadata.plau)
metadata


#data
##########
read_data <- function(data_filename, strain, wells_to_discard){
  data_raw <- read_excel(data_filename, skip=37)
  data_raw <- data_raw[-2,]
  
  data <- data_raw[-1, -1] %>% t() %>% as_tibble()
  names(data) <- data_raw[[1]][2:97]
  data <- data[-which(is.na(names(data)))]
  
  time <- data_raw[1,] %>% unlist()
  time <- time[-1] %>% as.numeric()
  data$time <- time/3600
  
  # limit the data until 36. hour (same time as other measurements)
  data <- data %>% filter(time <= 36)
  
  # discard contaminated/weird wells (identified earlier in separate files).
  data <- data %>%  dplyr::select(!all_of(wells_to_discard))
  
  # reshape data to long
  data_long <- 
    data %>% 
    pivot_longer(cols=A1:H12, names_to='well', values_to='OD')
  data_long$OD <- as.numeric(data_long$OD)
  
  #add column to identify strain
  data_long$strain <- strain
  
  return(data_long)
}

data.nema <- read_data('2025-01-14_Xnema/2025-01-14_nema.xlsx', 'X. nematophila', c("F12", "G12"))
data.bovi <- read_data('2025-02-06_Xbovi02/2025-02-06_Xbovi.xlsx', 'X. bovienii', c("G12"))
data.plau <- read_data('2025-01-31_Plaumo/2025-01-31_Plaumo.xlsx', 'P. laumondii', c("B1","C1","D1","E1","F1","G1","D12","E12","F12"))



data <- rbind(data.nema, data.bovi, data.plau)

# merge data & metadata

data_plot <- merge(data, metadata, by=c('strain','well'))
data_plot
#write_excel_csv(data_plot, file = 'all_data.csv')

# GROWTH CURVES
ggplot(data = data_plot %>% filter(strain!='blank'), 
       aes(x = time, y = OD, col=factor(concentration)))+
  #geom_point(size = 0.1)+
  geom_smooth(method = "loess", span = 0.5, se = F)+
  facet_wrap(~derivate+strain)+
  scale_colour_viridis_d()+
  ylim(c(0,0.5))+
  labs(x = "time [h]", col = "conc. [uM]")+
  theme_bw()
#theme(legend.position = "bottom")
  
  

####################################################################################### 

# Calculating the AUC

head(data_plot)

results <- 
  data_plot %>% 
  group_by(strain, well, precursor, derivate, concentration) %>%
  summarize(auc = auc(time, OD, ))

results$group <- paste(results$strain, results$derivate, sep='_')

ggplot(data = results,
       aes(x = concentration, y = auc, fill=factor(concentration), col=factor(concentration)))+
  geom_boxplot(width=0.5, staplewidth=0.5)+
  #geom_point()+
  facet_wrap(~derivate+strain)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  #xlim(c(0,24))+
  labs(x = "concentration [uM]", fill = "conc. [uM]")+
  guides(col='none')+
  #scale_x_log10()+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0,10,100,1000,2500))+
  theme_bw()
#theme(legend.position = "bottom")




### DRM 
groups <- unique(results$group)

fitlist <- list(groups)

for(group in groups){
  loopdata <- results[results$group == group,]
  cat("group", group, 
      "(", which(groups == group), "of", length(groups), "):")
  
  # pokud izolat neroste, preskocit
  # if(max(loopdata$OD, na.rm = TRUE) < 0.5){
  #   cat("NO GROWTH\n")
  #   
  #   next
  # }
  
  # ################### OUTLIERS
  # fit = drm(OD ~ koncentrace,
  #           data = loopdata,
  #           curveid = izolat,
  #           fct = W1.4(),
  #           na.action = na.omit,
  #           control = drmc(errorm = FALSE))
  # fit_bc <- boxcox(fit, plotit = FALSE) #nepot�ebuju BC, asi. Rezidua nemaj m�t norm�ln� rozd�len�, ne?
  # #EDIT: to asi ne, ale maj b�t homogenn�, nebo n�co. Homoskedastick�?
  # 
  # # �tverce odchylek maj� Chi2 rozd�len�
  # # spo��t� logick� vektor, kter� hodnoty jsou se zadanou prst� odlehl�
  # is_outlier <- scores(residuals(fit_bc)^2, type = "chisq", prob = 0.95)
  # 
  # filtered_data <- loopdata[!is_outlier,]
  
  fit = drm(auc ~ concentration,
              data = loopdata,
              curveid = group,
              #separate = TRUE,
              #pmodels = data.frame(izolat, 1, izolat, izolat), #rika, ze vsechny krivky maj stejny minimum
              fct = LL.4(),
              na.action = na.omit,
              control = drmc(errorm = FALSE))
  #fit_boxcox <- boxcox(fit, plotit = FALSE)
  
  
  fitlist[[group]] <- fit
  
  cat("DONE\n")
}

# IC50s
#takhle blb� to ned�l�m �pln� naschv�l, ale drc jako obvykle neni moc apply-friendly
groups <- groups %>% sort()

ED50list <- list()
for(group in groups){
  if(!is.null(fitlist[[group]])){
    ED50 <- ED(object = fitlist[[group]], 50)
    ED50 <- data.frame(ED50)
    ED50$group <- group
    
    ED50list <- rbind(ED50list, ED50)
  }
}
ED50list



# GRAF
# Let's only plot the ITCs, cause CNs are not really inhibited completely
groups <- groups[c(2,4,6)]
groups

par(mfrow = c(3,1),
    mar = c(2,3,0,0),
    oma = c(2,2,1,1))

for(group in groups){
  plot(fitlist[[group]], type = "all", legendPos = c(2, 2))+
    #text(x = 0.5, y = fitlist[[group]]$coefficients[3], labels = group)
    text(x = 1.5, y = 5, labels = group)
  abline(v = ED50list$Estimate[which(ED50list$group == group)], col = "red")
  print(group)
}

mtext("Compound concentration [uM]", side = 1, outer = TRUE)
mtext("auc", side = 2, outer = TRUE)




################################################################################
#
# Using normalized data
#
################################################################################

ggplot(data = results,
       aes(x = concentration, y = auc, fill=factor(concentration), col=factor(concentration)))+
  geom_boxplot(width=0.5, staplewidth=0.5)+
  #geom_point()+
  facet_wrap(~derivate+strain)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  #xlim(c(0,24))+
  labs(x = "concentration [uM]", fill = "conc. [uM]")+
  guides(col='none')+
  #scale_x_log10()+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0,10,100,1000,2500))+
  theme_bw()
#theme(legend.position = "bottom")

head(results)

# blank the results
results.bl <- 
  results %>% 
  group_by(strain, derivate, concentration) %>%
  mutate(auc_min = mean(auc) * (concentration==2500)) %>%
  group_by(strain, derivate) %>%
  mutate(auc_bl = auc - max(auc_min))

summary(results.bl)
results.bl %>% filter(derivate=='ITC', concentration==2500)

ggplot(data = results.bl,
       aes(x = concentration, y = auc_bl, fill=factor(concentration), col=factor(concentration)))+
  geom_boxplot(width=0.5, staplewidth=0.5)+
  #geom_point()+
  facet_wrap(~derivate+strain)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  #xlim(c(0,24))+
  labs(x = "concentration [uM]", fill = "conc. [uM]")+
  guides(col='none')+
  #scale_x_log10()+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0,10,100,1000))+
  theme_bw()
#theme(legend.position = "bottom")






# normalize the max of data to ~1
results.norm <- 
  results.bl %>% 
  group_by(strain, derivate, concentration) %>%
  mutate(auc_max = mean(auc_bl) * (concentration==0)) %>%
  group_by(strain, derivate) %>%
  mutate(auc_norm = auc_bl / max(auc_max))
  
  
summary(results.norm)
results.norm$auc_max
results.norm %>% filter(derivate=='ITC', concentration==0)



ggplot(data = results.norm %>% filter(derivate == 'ITC'),
       aes(x = concentration, y = auc_norm, fill=factor(concentration), col=factor(concentration)))+
  geom_boxplot(width=0.5, staplewidth=0.5)+
  #geom_point()+
  facet_wrap(~derivate+strain)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  #xlim(c(0,24))+
  labs(x = "concentration [uM]", fill = "conc. [uM]")+
  guides(col='none')+
  #scale_x_log10()+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0,10,100,1000))+
  theme_bw()
#theme(legend.position = "bottom")



# DRC
groups <- unique(results.norm$group)

fitlist <- list(groups)

for(group in groups){
  loopdata <- results.norm[results.norm$group == group,]
  cat("group", group, 
      "(", which(groups == group), "of", length(groups), "):")
  fit = drm(auc_norm ~ concentration,
            data = loopdata,
            curveid = group,
            #separate = TRUE,
            #pmodels = data.frame(izolat, 1, izolat, izolat), #rika, ze vsechny krivky maj stejny minimum
            fct = LL.4(),
            na.action = na.omit,
            control = drmc(errorm = FALSE))
  fitlist[[group]] <- fit
  
  cat("DONE\n")
}

# IC50s
#takhle blb� to ned�l�m �pln� naschv�l, ale drc jako obvykle neni moc apply-friendly
groups <- groups %>% sort()

ED50list <- list()
for(group in groups){
  if(!is.null(fitlist[[group]])){
    ED50 <- ED(object = fitlist[[group]], 50)
    ED50 <- data.frame(ED50)
    ED50$group <- group
    
    ED50list <- rbind(ED50list, ED50)
  }
}
ED50list



# GRAF
# Let's only plot the ITCs, cause CNs are not really inhibited completely
groups <- groups[c(2,4,6)]
groups

par(mfrow = c(3,1),
    mar = c(2,3,0,0),
    oma = c(2,2,1,1))

for(group in groups){
  plot(fitlist[[group]], type = "bar")+
    text(x = fitlist[[group]]$coefficients[4], y = fitlist[[group]]$coefficients[3], pos = 4, 
         labels = paste('IC50 =', 
                        sprintf(fitlist[[group]]$coefficients[4],fmt='%.1f'),
                                     'uM'))+
    text(x = 1, 
         y = (fitlist[[group]]$coefficients[3] + fitlist[[group]]$coefficients[2]) /2,
         pos=4, labels = group)
  abline(v = ED50list$Estimate[which(ED50list$group == group)], col = "red")
  print(group)
}

mtext("Compound concentration [uM]", side = 1, outer = TRUE)
mtext("AUC normalized [ ]", side = 2, outer = TRUE)
