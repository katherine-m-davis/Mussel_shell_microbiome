#### Mussel Microbiome Statistical Analyses ############################################################
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(forcats)
library(pairwiseAdonis)
library(indicspecies)
library(ggpubr)
beach.colors <- c("deepskyblue", "dodgerblue3", "goldenrod1", "sienna1", "grey55")
########################################################################################################################

#### Load data #####################################################################################################
# phyloseq
mmphy <- readRDS(file="~/Desktop/Musselmicrobe_phyloseq_filtered_rarefied_Feb2022.RDS")

# dataframe w/ updated chloroplast assignments
mmphydf <- readRDS(file="~/Desktop/Musselmicrobe_df_filtered_rarefied_March2022.RDS")
mmphydf$Tmt <- factor(mmphydf$Tmt)
mmphydf$Tmt <- fct_relevel(mmphydf$Tmt, "live control", "live", "killed")
########################################################################################################################


#### Summarise Shell erosion between paired treatments (add by Site) ###########################################################
mm.ae <- mmphydf %>% select(correct.name, Tmt, Site, area.er, Block.y) %>% distinct() %>% filter(!is.na(area.er))

sae <- ggpaired(mm.ae, x = "Tmt", y = "area.er", id = "Block.y", line.color = "gray", line.size = 0.4, fill = "Tmt") +
  stat_compare_means(label.x = 1.3, label.y = 80) +
  xlab('')+
  ylab('Shell area eroded (%)') +
  scale_fill_manual(values = c("grey22", "grey88")) +
  theme(legend.position = "none") #+
#facet_grid(~ Site)
sae
########################################################################################################################

#### Summarise cyanobacteria + chloroplast reads by paired treatments (add by Site) ############################################
mmphydf$Tmt <- factor(mmphydf$Tmt)
mmphydf$Tmt <- fct_relevel(mmphydf$Tmt, "live control", "live", "killed")
mmdf.co <- mmphydf %>% filter(Rank3 %in% c("Oxyphotobacteria", "Bacillariophyceae", "Coscinodiscophyceae","Mediophyceae", "Ulvophyceae", "unknown_Cyanobacteria")) %>% filter(!is.na(area.er))# Get cyanos
mmdf.co.sum <- mmdf.co %>% group_by(correct.name, Tmt, Site, area.er, Block.y) %>% summarise(c.abund = sum(Abundance)) %>% mutate(c.relabund = c.abund/1763 *100)

p <- ggpaired(mmdf.co.sum, x = "Tmt", y = "c.relabund", id = "Block.y", line.color = "gray", line.size = 0.4, fill = "Tmt") +
  stat_compare_means(label.x = 1.3, label.y = 80) +
  xlab('')+
  ylab('Relative abundance (%)') +
  scale_fill_manual(values = c("grey22", "grey88")) +
  theme(legend.position = "none", text = element_text(size=10)) +
  ylim(0,80)
########################################################################################################################

#### Summarise endolithic cyanobacteria by paired treatment (add by Site) ############################################
sort(unique(mmdf.co$Rank6))
pls <- c("unknown_Xenococcaceae","Pleurocapsa_PCC-7319", "Chroococcidiopsis_PCC-6712", "unknown_Nostocales")
mmdf.p <- mmphydf %>% filter(Rank6 %in% pls)
mmdf.p.sum <- mmdf.p %>% group_by(correct.name, Tmt, Site, elevation, area.er, Block.y, Rank6) %>% summarise(p.abund = sum(Abundance)) %>% mutate(p.relabund = p.abund/1763 *100)
mmdf.p.sum <- mmdf.p.sum %>% filter(! Tmt == "live control") %>% filter(! is.na(area.er))
mmdf.p.sum2 <- mmdf.p.sum %>% filter(p.relabund > 0)


e <- ggpaired(mmdf.p.sum2, x = "Tmt", y = "p.relabund", id = "Block.y", line.color = "gray", line.size = 0.4, fill = "Tmt") +
  stat_compare_means(label.x = 1.3, label.y = 80) +
  xlab('')+
  ylab('Relative abundance (%)') +
  scale_fill_manual(values = c("grey22", "grey88")) +
  theme(legend.position = "none", text = element_text(size=10)) +
  ylim(0,80)
e
#######################################################################################################################

#### Create combined plot of erosion and cyano relative abundance for paired treatments ################################
figure <- ggarrange(
  sae,                # First row with line plot
  # Second row with box and dot plots
  ggarrange(p, e, nrow = 2, labels = c("B", "C")), 
  nrow = 1, 
  labels = "A"       # Label of the line plot
) 
figure
ggsave("~/Desktop/paired_erosion_ra.png", width = 7, height = 6)
########################################################################################################################


#### Family level summaries for barplots ##############################################################################################
mmdf.fam <- mmphydf %>% group_by(Rank4, correct.name, Tmt, Site, elevation) %>% summarise(abund= sum(Abundance)) 
mmdf.fam <- mmdf.fam %>% mutate(relabund = abund/1763 * 100) # relative abundance

#Group by Beach and treatment and get families that are > 1.5 % abundant
mmdf.fam.sum  <-mmdf.fam %>% ungroup() %>% group_by(Site, Tmt, Rank4) %>% summarise(btabund = sum(abund))
mmdf.fam.sum <-mmdf.fam.sum %>% ungroup() %>% group_by(Site, Tmt) %>% mutate(totabund = sum(btabund))
mmdf.fam.sum <- mmdf.fam.sum %>% mutate(relabund = btabund/totabund * 100)

# Get Families that are > 1.5% of total per treatment + beach
mmdf.fam.sum.filt <- mmdf.fam.sum %>% filter(relabund > 1.5)
sort(unique(mmdf.fam.sum.filt$Rank4))

# Convert to "Other" if not present in top % families
mmdf.fam2 <- mmdf.fam
mmdf.fam2$Rank4o <- ifelse(mmdf.fam2$Rank4 %in% unique(mmdf.fam.sum.filt$Rank4), mmdf.fam2$Rank4, "Other")
unique(mmdf.fam2$Rank4o)

mmdf.fam2$Rank4o <- ifelse(mmdf.fam2$Rank4o %in% c("Eupodiscales", "Fragilariales", "unknown_green_algal_chloroplast"), "Chloroplast", mmdf.fam2$Rank4o)
unique(mmdf.fam2$Rank4o)


# Order samples by elevation
mmdf.fam2 <- mmdf.fam2 %>%
  mutate(correct.name = fct_reorder(correct.name, elevation))  
mmdf.fam2$elevation


# Reorder taxa 
dput(unique(mmdf.fam2$Rank4o))
mmdf.fam2$Rank4o <- factor(mmdf.fam2$Rank4o)
mmdf.fam2$Rank4o <- fct_relevel(mmdf.fam2$Rank4o, "Other","Chloroplast", "Nostocales", "Phormidesmiales", "Ardenticatenales", "Caulobacterales", "Chitinophagales", 
                                "Chromatiales", "Cytophagales", "Flavobacteriales", 
                                "Microtrichales", "Oceanospirillales", 
                                "Pirellulales", "Rhizobiales", "Rhodobacterales", "Rhodothermales", 
                                "Sphingomonadales", "Thiohalorhabdales", "Thiotrichales", "Tistrellales", 
                                "unknown_Gammaproteobacteria")
unique(mmdf.fam2$Rank4o)

# Reorder Treatments
mmdf.fam2$Tmt <- factor(mmdf.fam2$Tmt)
mmdf.fam2$Tmt <- fct_relevel(mmdf.fam2$Tmt, "live control", "live", "killed")

# Read in color paletter corresponding to Families
fc <- read.csv(file = "~/Desktop/MusselMicrobes/sequencecode/alyssa_mussel_dataset/KD_outputs/KD_data/Family_colors.csv")
fc <- read.csv(file = "~/Desktop/Family_colors.csv")
scales::show_col(fc$f.colors)
########################################################################################################################

#### Make barplots of dominant Families #################################################################################
ggplot(data=mmdf.fam2, aes(x=correct.name, y=relabund, fill=Rank4o)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  theme_classic() +
  facet_wrap(Site~Tmt, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = fc$f.colors) +
  theme(legend.title = element_blank(), legend.text=element_text(size=8)) +
  theme(legend.key.size = unit(0.3, 'cm'), text = element_text(size = 8)) +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") +
  ylab("Relative Abundance (%)") +
  guides(fill=guide_legend(ncol=1))
ggsave("~/Desktop/mm_barplots_family.png", width = 10, height = 8)
########################################################################################################################

#### Summarise cyanobacteria with new taxonomy assignments ##################################################
mmdf.c <- mmdf %>% group_by(Rank3,Rank4, Rank5, Rank6, correct.name, Tmt, Site, elevation, area.er, sequence) %>% summarise(abund= sum(Abundance)) 
mmdf.c <- mmdf.c %>% mutate(relabund = abund/1763 * 100) # relative abundance
# Filter taxa
mmdf.co <- mmdf.c %>% filter(Rank3 %in% c("Oxyphotobacteria", "Bacillariophyceae", "Coscinodiscophyceae","Mediophyceae", "Ulvophyceae", "unknown_Cyanobacteria")) # Get cyanos
sort(unique(mmdf.co$Rank5))
mmdf.co <- mmdf.co %>% ungroup()
mmdf.co$ASV <- paste("ASV", as.numeric(factor(mmdf.co$sequence)), sep = "_")
mmdf.co$ASV 
mmdf.co$Rank5.a <-paste(mmdf.co$Rank5, mmdf.co$ASV, sep = " ")

mmdf.co$Rank5.an <- ifelse(mmdf.co$Rank5 == "unknown_Chloroplast", mmdf.co$Rank5.a, mmdf.co$Rank5)
unique(mmdf.co$Rank5.an)

#Group by Beach and get taxa that are > 0.5 % abundant
mmdf.c.sum  <- mmdf.co %>% ungroup() %>% group_by(Site, Tmt, Rank5.an) %>% summarise(btabund = sum(abund))
mmdf.c.sum <-mmdf.c.sum %>% ungroup() %>%  group_by(Site, Tmt)%>%  mutate(totabund = sum(btabund))
mmdf.c.sum <- mmdf.c.sum %>% mutate(relabund = btabund/totabund * 100)

mmdf.c.sum.filt <- mmdf.c.sum %>% filter(relabund > 1.5)
dput(sort(unique(mmdf.c.sum.filt$Rank5.an)))

ucs <- c("unknown_Chloroplast ASV_126", "unknown_Chloroplast ASV_135", 
         "unknown_Chloroplast ASV_140", "unknown_Chloroplast ASV_289", 
         "unknown_Chloroplast ASV_300")
ucs.seq <- mmdf.co %>% ungroup() %>% filter(Rank5.an %in% ucs) %>% select(Rank5.an, sequence) %>% distinct()

ucs.seq$Rank6 <- NA
ucs.tax <- c("Haslea","Haslea","unknown Odontellaceae","unknown Bacillariophyta", "Polysiphonia")
ucs.seq$Rank6 <- ucs.tax

# Replace taxa assignments in the main dataframe
mmdf.co$Rank5.an <- ifelse(mmdf.co$Rank5.an %in% ucs.seq$Rank5.an, ucs.seq$Rank6, mmdf.co$Rank5.an)
unique(mmdf.co$Rank5.an)

# Replace taxa assignments in the subset dataframe
mmdf.c.sum.filt$Rank5.an <- ifelse(mmdf.c.sum.filt$Rank5.an %in% ucs.seq$Rank5.an, ucs.seq$Rank6, mmdf.c.sum.filt$Rank5.an)
unique(mmdf.c.sum.filt$Rank5.an)
########################################################################################################################

#### Filter to abundant cyanobacteria + chloroplast sequences ##################################################
# Convert to "Other" if not present in top % families
mmdf.co2 <- mmdf.co
mmdf.co2$Rank5o <- ifelse(mmdf.co2$Rank5.an %in% unique(mmdf.c.sum.filt$Rank5.an), mmdf.co2$Rank5.an, "Other")
mmdf.co2$Rank6o <- ifelse(mmdf.co2$Rank5.an %in% unique(mmdf.c.sum.filt$Rank5.an), mmdf.co2$Rank6, "Other")
mmdf.56 <- mmdf.co2 %>% ungroup() %>%  select(Rank5, Rank5o, Rank6, Rank6o) %>% distinct()
sort(unique(mmdf.co2$Rank6o))
# Change Rank 6 assignment
mmdf.co2$Rank6o <- ifelse(mmdf.co2$Rank6o=="unknown_Chloroplast", mmdf.co2$Rank5o, mmdf.co2$Rank6o)

# Add family assignments also
mmdf.co2$Rank6o_l <-ifelse(! mmdf.co2$Rank6o == "Other", paste0("(",mmdf.co2$Rank4,")"," ",mmdf.co2$Rank6o), "Other")
sort(unique(mmdf.co2$Rank6o_l)) # 15
mmdf.co2$Rank6o_l <-ifelse(mmdf.co2$Rank6o_l == "(unknown_green_algal_chloroplast) unknown_green_algal_chloroplast", "(Chloroplast) unknown_Green_alga", mmdf.co2$Rank6o_l)
mmdf.co2$Rank6o_l <-ifelse(mmdf.co2$Rank6o_l == "(Chloroplast) unknown Bacillariophyta", "(Chloroplast) unknown_Bacillariophyta", mmdf.co2$Rank6o_l)
mmdf.co2$Rank6o_l <- str_replace(mmdf.co2$Rank6o_l, "[(]Eupodiscales[)] ", "(Chloroplast) ")
mmdf.co2$Rank6o_l <- str_replace(mmdf.co2$Rank6o_l, "[(]Melosirales[)] ", "(Chloroplast) ")
sort(unique(mmdf.co2$Rank6o_l))

# Filter more to top 15 most abundant
mmdf.co2.sum <- mmdf.co2 %>% ungroup() %>% group_by(Rank6o_l) %>% summarise(t.abund =sum(abund))
mmdf.co2.sum <-mmdf.co2.sum %>% filter(t.abund >300)
mmdf.co2$Rank6o_l <- ifelse(mmdf.co2$Rank6o_l %in% mmdf.co2.sum$Rank6o_l, mmdf.co2.sum$Rank6o_l, "Other")

# Reorder Taxa assignments 
dput(unique(mmdf.co2$Rank6o_l))
mmdf.co2$Rank6o_l <- factor(mmdf.co2$Rank6o_l)
levels(mmdf.co2$Rank6o_l)
mmdf.co2$Rank6o_l <- fct_relevel(mmdf.co2$Rank6o_l, "Other", "(Chloroplast) Asterionellopsis_glacialis", 
                                 "(Chloroplast) Costaria_costata", "(Chloroplast) Ectocarpus_siliculosus", 
                                 "(Chloroplast) Melosira", "(Chloroplast) Odontella", "(Chloroplast) Paralemanea_annulata", 
                                 "(Chloroplast) Planoglabratella_opercularis", "(Chloroplast) Porphyra_umbilicalis_(laver)", 
                                 "(Chloroplast) unknown_Green_alga", "(Nostocales) Chroococcidiopsis_PCC-6712", 
                                 "(Nostocales) Pleurocapsa_PCC-7319", "(Nostocales) unknown_Oscillatoriaceae", 
                                 "(Phormidesmiales) Phormidesmis_ANT.LACV5.1")



# Get colors for plotting
# c.fam <- levels(mmdf.co2$Rank6o)
# write.csv(c.fam, file = "~/Desktop/MusselMicrobes/sequencecode/alyssa_mussel_dataset/KD_outputs/KD_data/Cyano_Chloro_colors_new.csv")
c.fam.colors <- read.csv(file = "~/Desktop/MusselMicrobes/sequencecode/alyssa_mussel_dataset/KD_outputs/KD_data/Cyano_Chloro_colors_new_v2.csv")
scales::show_col(c.fam.colors$colors)
########################################################################################################################

#### Make barplot of cyanobacteria + chloroplast sequences with new taxonomy assignments ################################
# summarise abundance at Rank6o_l
mmdf.co2.ab <- mmdf.co2 %>% ungroup() %>% group_by(correct.name, Site, Tmt, Rank6o_l, elevation) %>% summarise(abund2 = sum(abund)) %>% mutate(relabund2 = abund2/1763 *100)

op <- ggplot(data=mmdf.co2.ab, aes(x=correct.name, y=relabund2, fill=Rank6o_l)) + 
  geom_bar(color = "black", stat="identity", position="stack") +
  theme_classic() +
  facet_grid( ~ Tmt + Site, scales = "free_x", space = "free") +
  scale_fill_manual(values = c.fam.colors$colors) +
  theme(legend.title = element_text( size=8), legend.text=element_text(size=8)) +
  theme(legend.key.size = unit(0.3, 'cm'), text = element_text(size = 8)) +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "") +
  aes(x = fct_reorder(correct.name, elevation)) + # low to high 
  theme(panel.spacing = unit(.1, "lines")) +
  # panel.border = element_rect(color = "black", fill = NA, size = 1), 
  # strip.background = element_rect(color = "black", size = 1)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  theme(legend.title = element_blank()) wu
op
ggsave(file="~/Desktop/MusselMicrobes/MusselMicrobes_Figures/Cyanobacteria_barplots.svg", plot=op, width=17, height=4)
########################################################################################################################

#### Envfit analysis of transplant samples #############################################################################
mmphy.meta <- as_tibble(unclass(sample_data(mmphy)))
# Select envfit variables
mmphy.meta.env <- mmphy.meta %>% select(c("correct.name", "Site", "Tmt", "elevation", "air90th", "water90th","area.er"))
mmphy.meta.env <- mmphy.meta.env %>% filter(! is.na(elevation)) %>% filter(! is.na(area.er)) %>% filter(! is.na(air90th))
# subset phyloseq object to filtered dataset
mmphy.env <- subset_samples(mmphy, correct.name %in% mmphy.meta.env$correct.name) # filter samples and otu table to ones with available metadata
# Create distance matrix
m.env.dist.bc <- distance(mmphy.env, method="bray", type="samples") #bray-curtis
# Ordinate
pcoa.m.env <- ordinate(mmphy.env, method= "PCoA", distance = "bray")
# Get overview plot
plot_ordination(mmphy.env, pcoa.m.env, type="samples", color = "Site", shape = "Tmt") + 
  theme_classic()  +
  geom_point(size =4) +
  scale_color_manual(values = beach.colors)

# Run envfit
m.fit <- envfit(pcoa.m.env$vectors, mmphy.meta.env, permutations = 999, na.rm = T)
m.fit
m.env.scores <- as.data.frame(scores(m.fit, display = "vectors"))
m.env.scores <- cbind(m.env.scores, pval = m.fit$vectors$pvals)
m.env.scores.sig <- m.env.scores %>% filter(pval < 0.05)
rownames(m.env.scores.sig) <- c("elevation", "subaerial T (90th quantile)", "immersed T (90th quantile)")
# 
#             Axis.1   Axis.2     r2 Pr(>r)    
# elevation -0.10055  0.99493 0.1434  0.028 *  
# air90th    0.24781 -0.96881 0.1298  0.036 *  
# water90th -0.95072  0.31007 0.3118  0.001 ***
# area.er    0.88550 -0.46464 0.0042  0.911  
# 
# Goodness of fit:
#               r2    Pr(>r)    
# correct.name 1.0000  1.000    
# Site         0.2361  0.001 ***
# Tmt          0.0173  0.457  
########################################################################################################################

#### Plot ordination with envfit explanatory arrows ####
# get PCoA coordinates
pcoa.axes <- as.data.frame(pcoa.m.env$vectors[,1:2])
pcoa.axes$correct.name <- rownames(pcoa.axes)
pcoa.axes <- left_join(pcoa.axes, mmphy.meta.env)  # add metadata
# Plot
ggplot(pcoa.axes) +  
  geom_point(mapping = aes(x = Axis.1, y = Axis.2, colour = Site, shape = Tmt), size =4) +
  xlab("Axis 1 [16.8%]") + 
  ylab("Axis 2 [11.3%]") +
  theme_classic() +
  coord_fixed() + ## need aspect ratio of 1
  scale_color_manual(values = beach.colors) +
  labs(colour = "Site", shape = "Treatment") +
  geom_segment(data= m.env.scores.sig, aes(x = 0, xend = Axis.1, y = 0, yend = Axis.2), arrow = arrow(length = unit(0.2, "cm")), color = "grey10", lwd =0.3) +
  ggrepel::geom_text_repel(data = m.env.scores.sig, aes(x = Axis.1, y= Axis.2+.045, label = rownames(m.env.scores.sig)))
ggsave("~/Desktop/mm_envfit_pcoa.png", width = 6, height = 5)
########################################################################################################################




#### Summarise Cyanobacteria reads by Treatment and Site ###################################################################################
mmdf.c <- mmphydf %>% group_by(Rank3,Rank4, Rank5, Rank6, correct.name, Tmt, Site, elevation, area.er) %>% summarise(abund= sum(Abundance)) 
mmdf.c<- mmdf.c %>% mutate(relabund = abund/1763 * 100) # relative abundance
mmdf.co <- mmdf.c %>% filter(Rank3 == "Oxyphotobacteria") # Get cyanos
# Remove Chloroplast sequences
unique(mmdf.co$Rank4)
mmdf.cnc <- mmdf.co %>% filter(! Rank4 == "Chloroplast")
mmdf.cnc.sum <- mmdf.cnc %>% ungroup() %>% group_by(correct.name, Tmt, Site, elevation, area.er) %>% summarise(cy.abund = sum(abund))
mmdf.cnc.sum <- mmdf.cnc.sum %>% ungroup() %>%  mutate(cy.relabund = cy.abund/1763 * 100)

