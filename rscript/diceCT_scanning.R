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

### Set working dir
### Type in  your own path in the quotation marks
dir <- "~/guide-diceCT-snakes/data/" # specimen list location
dir_gv <- "~/guide-diceCT-snakes/data/Raw_grayValue_data/" # grayvalue location
dir_plots <- "~/guide-diceCT-snakes/plots/"

### Load data 
# Loads in datasheet with specimen list (includes taxonomy, staining duration and size measurements)
dat <- read_csv(paste0(dir, "specimen_list.csv"), col_names = TRUE) %>%
  select(Taxon_family, Genus, Species, G_species, RAB = `RAB #`, Taxon_name, Taxon_ID, Museum, Specimen = `Specimen #`,
         SVL_mm, Mass_g, Head_diameter_mm = HeadGirth_mm, Days_stained = `Days Stained`) %>%
  mutate(Head_radius_mm = Head_diameter_mm/2, 
         Head_diffusion_rate = Days_stained/Head_radius_mm) %>%
  arrange(Taxon_family, Genus, Species)

### Summarise specimen data
dat %>% group_by(Taxon_family) %>% count() # total number of specimens per family
mean(dat$Head_diffusion_rate) # mean rate of iodine diffusion for the head
sd(dat$Head_diffusion_rate)
mean(dat$Body_diffusion_rate, na.rm = T) # mean rate of iodine diffusion for the body

### Relationship between staining and specimen size
## Linear regression
head <- lm(log(Days_stained) ~ log(Head_radius_mm), dat) # head radius linear model
mass <- lm(log(Days_stained) ~ log(Mass_g), dat) # mass
svl <- lm(log(Days_stained) ~ log(SVL_mm), dat) # svl

## Summary of lm
summary(head)
coef(head)  
confint(head)

## Plot staining duration over specimen size

# days~head
phead <- dat %>%
  ggplot(aes(x = Head_radius_mm, y = Days_stained)) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(method = "lm", alpha = .2, colour = "grey20") + # add linear model line
  scale_y_continuous(breaks = seq(min(dat$Days_stained), max(dat$Days_stained), by = 3), limits = c(3,12)) +
  scale_x_continuous(breaks = round(seq(min(dat$Head_radius_mm), max(dat$Head_radius_mm), by = 3), 0)) +
  labs(x = "Head radius (mm)", y = "") + # plot labels
  theme_half_open() +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.line = element_line(size = 1))

# days~mass
pmass  <- dat %>%
  ggplot(aes(x = log(Mass_g), y = Days_stained)) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(method = "lm", alpha = .2, colour = "grey20") + 
  scale_y_continuous(breaks = seq(min(dat$Days_stained), max(dat$Days_stained), by = 3), limits = c(3,12)) +
  labs(x = "ln Mass", y = "") + 
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
  labs(x = "SVL (mm)", y = "Staining duration (d)") + 
  theme_half_open() +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.line = element_line(size = 1))

# Combine plots to make complete Figure
all_plots <- plot_grid(psvl, pmass, phead, ncol = 3, labels = c("(a)", "(b)", "(c)"))
all_plots

# Save plot
ggsave(filename = "fig_5.png",
       device = "png",
       path = normalizePath(dir_plots),
       plot = all_plots,
       width = 297, 
       height = 87, 
       units = "mm")

### Plot grayscale values for diceCT scans as boxplot and histograms

## Load in raw gray value data
# I had to manually change some of the file names to make sure this worked
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
gv <- gv %>%
  full_join(dat)

## boxplot
pbox <- gv %>%
  filter(Voxel_count < 2500) %>%
  ggplot(aes(reorder(G_species, Voxel_count), Voxel_count, fill = Days_stained)) +
  geom_boxplot(shape = 1) +
  scale_fill_viridis_c(alpha=0.6, breaks = seq(3,12, by = 2), direction = -1, option = "inferno") +
  theme_classic() +
  labs(x = "", y = "", fill = "Days") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic", size = 25),
        axis.line.x = element_line(colour = 'black', size = 2),
        axis.ticks.x = element_line(colour = "black", size = 2),
        axis.line.y = element_line(colour = 'black', size = 2),
        axis.ticks.y = element_line(colour = "black", size = 2),
        text = element_text(size = 20))

pbox

# save boxplot
ggsave(filename = "boxplot.png",
       device = "png",
       path = normalizePath(dir_plots),
       plot = pbox,
       width = 520, 
       height = 150, 
       units = "mm")

## histograms
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
