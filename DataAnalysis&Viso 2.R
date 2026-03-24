# Meditation exercise: Endpoint-focused vs. Present-focused
# Behavioral data analysis & plotting
# For: emotional experiences
# Requires: "rawdata.xlsx" (sheet = "rawdata1")
# Programmed by Feng XIAO (updated on 2026-3-18)

####################################################################################################
### 0) Preparation ---------------------------------------------------------------------------------
####################################################################################################
# Load packages actually used
pkg <- c(
  "readxl","dplyr","tidyr","stringr","forcats",
  "ggplot2","janitor","car","withr","emmeans",
  "lme4","lmerTest","nnet"
)
lapply(pkg, require, character.only = TRUE)

# Working directory (RStudio only)
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Global ggplot defaults
theme_set(theme_classic(base_size = 8))

####################################################################################################
### 1) Data input & cleaning -------------------------------------------------------------
####################################################################################################

rd <- readxl::read_excel("rawdata.xlsx", sheet = "rawdata1", na = c("", "NA", "---")) %>%
  mutate(
    Group  = factor(Group, levels = c("present","endpoint")), # set present-focused as baseline
    Gender = factor(Gender, levels = c(1,2), labels = c("male","female")),
    Age_z  = as.numeric(scale(Age))
  )

emo_long <- rd %>%
  select(SubjNum, Group, Gender, Age_z,
         Emotion1, Emotion2, Emotion3,
         Rating1,  Rating2,  Rating3) %>%
  pivot_longer(
    cols = c(Emotion1, Emotion2, Emotion3, Rating1, Rating2, Rating3),
    names_to = c(".value","slot"),
    names_pattern = "(Emotion|Rating)([1-3])"
  ) %>%
  mutate(
    Emotion = Emotion %>%
      as.character() %>%
      str_trim() %>%
      na_if("") %>%
      str_to_title()
  ) %>%
  # delete NA
  filter(!is.na(Emotion), !is.na(Rating)) %>%
  # set present-focused as baseline
  mutate(
    Group = relevel(Group, ref = "present")
  )

# Check for emotion type and confirn the Top-K
# Criteria: report frequency over 30
emo_freq <- emo_long %>% count(Emotion, sort = TRUE)
print(emo_freq, n = 30) # Three emotions: Peace, Relaxation, Sadness

emo_table <- emo_long %>%
  count(Group, Emotion) %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(Total = present + endpoint) %>%
  arrange(desc(Total))
print(emo_table, n = nrow(emo_table))
readr::write_csv(emo_table, "Table1_EmotionType_byGroup_simple.csv")

## Focus on the 3 frequent emotions
K <- 3
top_emos <- emo_freq %>% slice_head(n = K) %>% pull(Emotion)
emo_long <- emo_long %>%
  mutate(Emotion_collapsed = if_else(Emotion %in% top_emos, Emotion, "Other"),
         Emotion_collapsed = fct_infreq(Emotion_collapsed)) 

####################################################################################################
### 2) Analysis for emotional types -------------------------------------------------------------
####################################################################################################

# lock order for readability; drop "Other"
emo3_levels <- c("Peace","Relaxation","Sadness")
emo3 <- emo_long %>%
  filter(Emotion %in% emo3_levels) %>%
  mutate(
    Emotion = factor(Emotion, levels = emo3_levels),
    # keep a clean variable name used below
    Emotion3 = Emotion
  )

# sanity check: counts by group ¡Á emotion
emo3 %>% count(Group, Emotion3) %>% tidyr::pivot_wider(names_from = Emotion3, values_from = n, values_fill = 0)

# Group difference in emotion-type distribution
# Chi-square test on the 3-class distribution
tab_GE3 <- table(emo3$Group, emo3$Emotion3)
chisq_res3 <- chisq.test(tab_GE3)
chisq_res3                     # overall test
chisq_res3$stdres              # standardized residuals (which cells drive the effect)

# Multinomial logistic regression:
# Emotion3 ~ Group + Age_z + Gender
# set the baseline as the most frequent emotion (here: "Peace")
emo3$Emotion3 <- relevel(emo3$Emotion3, ref = "Peace")
fit_multi3 <- nnet::multinom(Emotion3 ~ Group + Age_z + Gender, data = emo3, trace = FALSE)
summary(fit_multi3)

# Wald z-tests p-values for coefficients
zvals3 <- summary(fit_multi3)$coefficients / summary(fit_multi3)$standard.errors
pvals3 <- 2 * (1 - pnorm(abs(zvals3)))
pvals3

####################################################################################################
### 3) Analysis for emotional ratings -------------------------------------------------------------
####################################################################################################

## Ratings: Group ¡Á Emotion mixed model
# LMM: Rating ~ Group * Emotion3 + Age_z + Gender + (1|SubjNum)
m_rate3 <- lmer(Rating ~ Group * Emotion3 + Age_z + Gender + (1|SubjNum), data = emo3)

# Type-III tests (focus on Emotion3 and Group:Emotion3)
withr::with_options(
  list(contrasts = c("contr.sum","contr.poly")),
  print(car::Anova(m_rate3, type = 3))
)
coef(summary(m_rate3))[c("Age_z","Genderfemale"), ]

# Estimated marginal means and pairwise tests (Holm correction)
# (B1) Between-group at each emotion
emm_bt3 <- emmeans(m_rate3, ~ Group | Emotion3)
bt3_tbl <- summary(pairs(emm_bt3, adjust = "holm"), infer = c(TRUE, TRUE))
bt3_tbl

# (B2) Within-group across emotions
emm_wi3 <- emmeans(m_rate3, ~ Emotion3 | Group)
wi3_tbl <- summary(pairs(emm_wi3, adjust = "holm"), infer = c(TRUE, TRUE))
wi3_tbl

####################################################################################################
### 4) Plotting for emotional types and ratings -------------------------------------------------------------
####################################################################################################

# Emotional proportions
emo_prop <- emo3 %>%
  dplyr::count(Group, Emotion3, name = "k") %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(N = sum(k),
                prop = k / N) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(ci = list(prop.test(k, N, correct = FALSE)$conf.int),
                p_low = ci[1], p_high = ci[2]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Group   = factor(Group, levels = c("endpoint","present")),
    Emotion3 = factor(Emotion3, levels = c("Peace","Relaxation","Sadness"))
  )

# Aesthetics setting
grp_cols   <- c("present" = "#4169E1", "endpoint" = "#B22222")
pd_bar     <- position_dodge(width = 0.60)

p_emoprop <- ggplot(emo_prop, aes(x = Emotion3, y = prop, fill = Group)) +
  geom_col(position = pd_bar, width = 0.50,
           colour = "black", linewidth = 0.35, alpha = 0.70) +
  geom_errorbar(aes(ymin = p_low, ymax = p_high),
                position = pd_bar, width = 0.12, size = 0.45, colour = "black") +
  scale_fill_manual(values = grp_cols) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion within group (%)") +
  theme(
    axis.line   = element_line(colour = "black", linewidth = 0.35),
    axis.title  = element_text(size = 7, colour = "black"),
    axis.text   = element_text(size = 7, colour = "black"),
    legend.position = "none",
    plot.margin = ggplot2::margin(3, 3, 2, 2)
  ) +
  patchwork::plot_annotation(
    title = "(d) Emotion type distribution (top-3)",
    theme = theme(plot.title = element_text(size = 7, colour = "black", face = "bold",
                                            margin = ggplot2::margin(b = 2)))
  )

ggsave("pic_emoProp.pdf", plot = p_emoprop, width = 2, height = 2, units = "in",
       device = cairo_pdf)

####################################################################################################
# END
####################################################################################################
