### Guide to diceCT scanning snakes ###

### Script to 1) plot regression iodine staining duration / specimen size, and
### 2) plot grayscale values from diceCT scans

### Load packages
# Install packages before loading them using this function 'install.packages("packageName")'
library(tidyverse)
library(magrittr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(grid)
library(viridisLite)

### Set working dir
### Type in  your own path in the quotation marks
dir <- "C:/Users/jmcr/Documents/guide-diceCT-snakes/data/" # specimen list location
dir_gv <- "C:/Users/jmcr/Documents/guide-diceCT-snakes/data/Raw_grayValue_data/" # grayvalue location
dir_plots <- "C:/Users/jmcr/Documents/guide-diceCT-snakes/plots/"

### Load data 
# Loads in datasheet with specimen list (includes taxonomy, staining duration and size measurements)
dat <- read_csv(paste0(dir, "specimen_list.csv"), col_names = TRUE) %>%
  select(Taxon_family, Genus, Species, G_species, RAB = `RAB #`, Taxon_name, Taxon_ID, Museum, Specimen_ID = `Specimen #`,
         SVL_mm, Mass_g, Head_diameter_mm = HeadGirth_mm, Days_stained = `Days Stained`, Preservation_age) %>%
  mutate(Head_radius_mm = Head_diameter_mm/2, 
         Head_diffusion_rate = Days_stained/Head_radius_mm) %>%
  arrange(Taxon_family, Genus, Species)

### Summarise specimen data
dat %>% group_by(Taxon_family) %>% count() # total number of specimens per family
mean(dat$Head_diffusion_rate) # mean rate of iodine diffusion for the head
sd(dat$Head_diffusion_rate)

### Relationship between staining and specimen size
## Linear regression
head <- lm(log(Days_stained) ~ log(Head_radius_mm), dat) # head radius linear model
head_r2 <- round(summary(head)$r.squared, 2) # save r squared value

mass <- lm(log(Days_stained) ~ log(Mass_g), dat) # mass
mass_r2 <- round(summary(mass)$r.squared, 2)

svl <- lm(log(Days_stained) ~ log(SVL_mm), dat) # svl
svl_r2 <- round(summary(svl)$r.squared, 2)

age <- lm(log(Days_stained) ~ log(Preservation_age), dat)

## Summary of lm
summary(head)
coef(head)  
confint(head)

summary(age)
coef(age)
confint(age)

## Plot staining duration over specimen size

# days~head
phead <- dat %>%
  ggplot(aes(x = Head_radius_mm, y = Days_stained)) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(method = "lm", alpha = .2, colour = "grey20") + # add linear model line
  scale_y_continuous(breaks = seq(min(dat$Days_stained), max(dat$Days_stained), by = 3), limits = c(3,12)) +
  scale_x_continuous(breaks = round(seq(min(dat$Head_radius_mm), max(dat$Head_radius_mm), by = 3), 0)) +
  annotate("text", x = 2.5, y = 10, label = paste("italic(R) ^ 2 == ", head_r2), parse = TRUE) +
  labs(x = "Head radius (mm)", y = "Staining duration (d)") + # plot labels
  theme_half_open() +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.line = element_line(size = 1))

# days~mass
pmass  <- dat %>%
  ggplot(aes(x = log(Mass_g), y = Days_stained)) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(method = "lm", alpha = .2, colour = "grey20") + 
  scale_y_continuous(breaks = seq(min(dat$Days_stained), max(dat$Days_stained), by = 3), limits = c(3,12)) +
  annotate("text", x = 2, y = 10, label = paste("italic(R) ^ 2 == ", mass_r2), parse = TRUE) +
  labs(x = "ln Mass", y = "Staining duration (d)") + 
  theme_half_open() +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.line = element_line(size = 1))

# days~svl
psvl  <- dat %>%
  drop_na(SVL_mm) %>%
  ggplot(aes(x = SVL_mm, y = Days_stained)) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(method = "lm", alpha = .2, colour = "grey20") + 
  scale_y_continuous(breaks = seq(min(dat$Days_stained), max(dat$Days_stained), by = 3), limits = c(3,12)) +
  xlim(230, 1840) +
  annotate("text", x = 320, y = 10, label = paste("italic(R) ^ 2 == ", svl_r2), parse = TRUE) +
  labs(x = "SVL (mm)", y = "Staining duration (d)") + 
  theme_half_open() +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.line = element_line(size = 1))

# days~age
dat %>%
  drop_na(Preservation_age) %>%
  ggplot(aes(x = log(Preservation_age), y = Days_stained)) +
  geom_point(shape = 1, size = 2) +
  #geom_smooth(method = "lm", alpha = .2, colour = "grey20") + 
  scale_y_continuous(breaks = seq(min(dat$Days_stained), max(dat$Days_stained), by = 3), limits = c(3,12)) +
  #scale_x_continuous(breaks = seq(min(0.5), max(5), by = 0.5), limits = c(0.5,5)) +
  #xlim(230, 1840) +
  labs(x = "Specimen age (years)", y = "Staining duration (d)") + 
  theme_half_open() +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.line = element_line(size = 1))

# Combine plots to make complete Figure
all_plots <- plot_grid(phead, psvl, pmass, ncol = 1, labels = c("(a)", "(b)", "(c)"))
all_plots

# Save plot
ggsave(filename = "fig_3.png",
       device = "png",
       path = normalizePath(dir_plots),
       plot = all_plots,
       width = 125, 
       height = 295, 
       units = "mm")

### Plot grayscale values for diceCT scans as boxplot and histograms

## Load in raw gray value data
gv <- list.files(path = dir_gv, pattern = ".csv", full.names = TRUE) %>%
  set_names(str_remove(string = basename(.), pattern = "raw.csv")) %>%
  map(read.csv, stringsAsFactors = FALSE) %>%
  bind_rows(.id = "sample") %>%
  separate(col = "sample", into = c("Museum", "Collection", "Specimen_ID"), sep = "-", extra = "merge") %>%
  separate(col = "Specimen_ID", into = c("Specimen_ID", "Genus_species"), sep = "_", extra = "drop") %>%
  #separate(col = "Genus_species", into = c("Genus", "Species"), sep = "-", extra = "merge") %>%
  select(Museum, Specimen_ID, Genus_species, Gray_value = `ï..Gray.value`, Voxel_count = `Voxel.count`) %>%
  mutate(Gray_value = str_remove(string = Gray_value, pattern = "> ")) %>%
  mutate(Genus_species = str_replace(string = Genus_species, pattern = "-", replacement = "_")) %>%
  mutate(Museum = str_to_upper(Museum))

# Change char class to numeric for some columns
gv$Gray_value <- as.numeric(gv$Gray_value)
gv$Gray_value <- round(gv$Gray_value)
gv$Specimen_ID <- as.numeric(gv$Specimen_ID)

# Join gray value data with specimen info
box_p <- subset(gv, Voxel_count <2500) 
box_p <- box_p %>%
  full_join(dat) 

head(box_p)
tail(box_p)

## boxplot
bp <- box_p %>%
  ggplot(aes(reorder(Taxon_ID, Voxel_count), Voxel_count, fill = Days_stained)) +
  geom_boxplot(shape = 1) +
  scale_fill_viridis_c(alpha=0.6, breaks = seq(3,12, by = 2), direction = -1, option = "inferno") +
  #geom_text(aes(label = Days_stained), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_classic() +
  labs(x = "", y = "", fill = "Days") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic", size = 25),
        axis.line.x = element_line(colour = 'black', size = 2),
        axis.ticks.x = element_line(colour = "black", size = 2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.ticks.y = element_line(colour = "black", size = 2),
        text = element_text(size = 20))


# save boxplot
ggsave(filename = "boxplot2.png",
       device = "png",
       path = normalizePath(dir_plots),
       plot = bp,
       width = 520, 
       height = 150, 
       units = "mm")

## histograms
# Change name of species to plot different species
hist <- gv %>%
  filter(Genus_species %in% "Bothrops_bilineatus", Voxel_count < 100000) %>%
  ggplot(aes(x = Gray_value, y = Voxel_count)) +
  geom_bar(stat = "identity", width = 1, alpha = 0.5) +
  theme_classic() +
  scale_y_continuous() +
  labs(x = "Gray value", y = "Voxel count", title = "") +
  theme(axis.line.x = element_line(colour = 'black', size = 2),
        axis.ticks.x = element_line(colour = "black", size = 2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.ticks.y = element_line(colour = "black", size = 2),
        text = element_text(size = 40))

hist

# save boxplot
ggsave(filename = "Bobi_UMMZ245084.png",
       device = "png",
       path = normalizePath(dir_plots),
       plot = hist,
       width = 297, 
       height = 210,
       units = "mm")

###
sessionInfo()
###
