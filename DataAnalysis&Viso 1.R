# Mental exercise: Endpoint-focused vs. Present-focused
# Behavioral data analysis & plotting
# For: delay discounting, engagement, time perception, resource allocation
# Requires: "rawdata.xlsx" (sheet = "rawdata1")
# Programmed by Feng XIAO (updated on 2025-10-25)

####################################################################################################
### 0) Preparation ---------------------------------------------------------------------------------
####################################################################################################

# Load packages actually used
pkg <- c(
  "car","tidyr","dplyr","readxl","effectsize","lme4","lmerTest",
  "emmeans","ggplot2","patchwork","withr","randomForest"
)
lapply(pkg, require, character.only = TRUE)

# Working directory (RStudio only)
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Global ggplot defaults (journal-like)
theme_set(theme_classic(base_size = 8))

# Helper: partial eta-squared from car::Anova(Type III) table
eta_from_Anova <- function(aov_tbl, term = "Group") {
  DF1 <- as.numeric(aov_tbl[term, "Df"])
  Fv  <- as.numeric(aov_tbl[term, "F value"])
  DF2 <- as.numeric(aov_tbl["Residuals", "Df"])
  effectsize::F_to_eta2(Fv, df = DF1, df_error = DF2, partial = TRUE, ci = 0.95)
}

# Helper: identify Tukey outliers for a given numeric vector
flag_outliers <- function(x) {
  q1  <- quantile(x, 0.25, na.rm = TRUE)
  q3  <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  (x < q1 - 1.5*iqr) | (x > q3 + 1.5*iqr)
}

####################################################################################################
### 1) Data input & quick demographics -------------------------------------------------------------
####################################################################################################

rd <- readxl::read_excel("rawdata.xlsx", sheet = "rawdata1", na = "---")
rd_endpoint <- filter(rd, Group == 'endpoint')
rd_present <- filter(rd, Group == 'present')

# Quick demographics summary
demo <- rd %>%
  mutate(Gender = factor(Gender, levels = c(1,2), labels = c("male","female")),
         Group  = factor(Group,  levels = c("endpoint","present"))) %>%
  summarise(
    n           = n(),
    n_male      = sum(Gender == "male", na.rm = TRUE),
    n_female    = sum(Gender == "female", na.rm = TRUE),
    age_mean    = mean(Age, na.rm = TRUE),
    age_sd      = sd(Age, na.rm = TRUE)
  )

# Overall
# Age
mean(rd$Age) #26.61
sd(rd$Age) #8.49
mean((filter(rd, Gender == 1))$Age) #male: 28.96
sd((filter(rd, Gender == 1))$Age) #male: 8.11
mean((filter(rd, Gender == 2))$Age) #female: 25.38
sd((filter(rd, Gender == 2))$Age) #female: 8.47
t.test((filter(rd, Gender == 1))$Age, (filter(rd, Gender == 2))$Age,
       paired =FALSE, alternative = c("two.sided"), var.equal=FALSE,
       conf.level=0.95) #Age: Male>female, p=.014

# Gender
dim(filter(rd, Gender == 1))[1] #51 males
dim(filter(rd, Gender == 2))[1] #97 females

# Endpoint-focused group
# Age
mean(rd_endpoint$Age) #26.72
sd(rd_endpoint$Age) #8.69
mean((filter(rd_endpoint, Gender == 1))$Age) #male: 29.70
sd((filter(rd_endpoint, Gender == 1))$Age) #male: 9.16
mean((filter(rd_endpoint, Gender == 2))$Age) #female: 24.68
sd((filter(rd_endpoint, Gender == 2))$Age) #female: 7.82

# Gender
dim(filter(rd_endpoint, Gender == 1))[1] #30 males
dim(filter(rd_endpoint, Gender == 2))[1] #44 females

# Present-focused group
# Age
mean(rd_present$Age) #26.51
sd(rd_present$Age) #8.35
mean((filter(rd_present, Gender == 1))$Age) #male: 27.90
sd((filter(rd_present, Gender == 1))$Age) #male: 6.39
mean((filter(rd_present, Gender == 2))$Age) #female: 25.96
sd((filter(rd_present, Gender == 2))$Age) #female: 9.01

# Gender
dim(filter(rd_present, Gender == 1))[1] #21 males
dim(filter(rd_present, Gender == 2))[1] #53 females

####################################################################################################
### 2) LMM: Delay discounting  ---------------------------------------------------------------------
####################################################################################################

# Data preparation
rd <- rd %>%
  mutate(
    logk_pre   = log(k_pre),
    logk_post1 = log(k_post1),
    logk_post2 = log(k_post2),
    logk_post3 = log(k_post3),
    Gender = factor(Gender, levels = c(1,2), labels = c("male","female")),
    Group  = relevel(factor(Group, levels = c("endpoint","present")), ref = "present"),
    Age_z  = as.numeric(scale(Age))
  )

rd_long <- rd %>%
  pivot_longer(starts_with("logk_"), names_to = "time", values_to = "logk",
               names_prefix = "logk_") %>%
  filter(!is.na(logk)) %>%
  mutate(time = factor(time, levels = c("pre","post1","post2","post3")))

# Model (random intercept by subject)
model_logk <- lmer(logk ~ Group * time + Age_z + Gender + (1 | SubjNum),
                   data = rd_long, REML = TRUE)

# Type III tests (local contrasts setting)
withr::with_options(
  list(contrasts = c("contr.sum","contr.poly")),
  print(car::Anova(model_logk, type = 3))
)

# Pairwise (between groups at each time; Holm)
emm_bt <- emmeans(model_logk, ~ Group | time)
bt_tbl <- summary(contrast(emm_bt, "pairwise", adjust = "holm"), infer = c(TRUE, TRUE))
bt_tbl

# Pairwise (within group across times; Holm)
emm_wi <- emmeans(model_logk, ~ time | Group)
wi_tbl <- summary(contrast(emm_wi, "pairwise", adjust = "holm"), infer = c(TRUE, TRUE))
wi_tbl

# Effect size (Cohen's d): between-group at each time (residual SD)
d_bt <- eff_size(emm_bt, sigma = sigma(model_logk), edf = df.residual(model_logk)) %>%
  as.data.frame() %>%
  select(time, contrast, d = effect.size, SE, df, lower.CL, upper.CL)
d_bt

# Effect size (Cohen's d): within-group pairwise changes (residual SD)
chg_tbl <- contrast(emm_wi, "pairwise") %>% summary(infer = c(TRUE, TRUE)) %>%
  as.data.frame() %>%
  mutate(d = estimate / sigma(model_logk))
chg_tbl

# Plot: adjusted means over time
em_plot <- emmeans(model_logk, ~ time * Group) %>%
  as.data.frame() %>%
  mutate(Group = factor(Group, levels = c("endpoint","present")),
         time  = factor(time,  levels = c("pre","post1","post2","post3")))

group_colors <- c("present" = "#4169E1", "endpoint" = "#B22222")
group_shapes <- c("endpoint" = 21, "present" = 24)
pd <- position_dodge(width = 0.30)

p_dd <- ggplot(em_plot, aes(time, emmean, group = Group, colour = Group, shape = Group)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.10, linewidth = 0.45, position = pd) +
  geom_point(position = pd, fill = "white", size = 1.0, stroke = 0.7) +
  geom_line(position = pd, linewidth = 0.45, linetype = "solid") +
  scale_colour_manual(values = group_colors) +
  scale_shape_manual(values  = group_shapes) +
  scale_y_continuous(limits = c(-8, -4), breaks = seq(-8, -4, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Log k") +
  theme(
    axis.line   = element_line(colour = "black", linewidth = 0.35),
    axis.title  = element_text(size = 7, colour = "black"),
    axis.text   = element_text(size = 7, colour = "black"),
    legend.position = "none",
    plot.margin = ggplot2::margin(3, 3, 2, 2)
  ) +
  patchwork::plot_annotation(
    title = "(a) Delay discounting",
    theme = theme(plot.title = element_text(size = 7, colour = "black", face = "bold",
                                            margin = ggplot2::margin(b = 2)))
  )

ggsave("pic_dd.pdf", plot = p_dd, width = 2, height = 2, units = "in", device = cairo_pdf)

####################################################################################################
### 3) ANCOVA: Subjective time perception ----------------------------------------------------------
####################################################################################################

dat <- rd %>%
  select(SubjNum, Group, Gender, Age, TenPercp,
         Account_Immediate, Account_ShortTerm, Account_LongTerm) %>%
  filter(!is.na(Age), !is.na(Gender), !is.na(Group)) %>%
  mutate(Gender = factor(Gender, levels = c("male","female")),
         Group  = relevel(factor(Group), ref = "present"),
         Age_z  = as.numeric(scale(Age)))

m_time <- lm(TenPercp ~ Group + Age_z + Gender, data = dat)

a3_time <- withr::with_options(
  list(contrasts = c("contr.sum","contr.poly")),
  car::Anova(m_time, type = 3)
)
print(a3_time)

eta_from_Anova(as.data.frame(a3_time), "Group")

# Interaction check
anova(m_time, lm(TenPercp ~ Group*Gender + Group*Age_z + Age_z + Gender, data = dat))

# Plot (same half-split style)
dat <- dat %>%
  mutate(Group = factor(Group, levels = c("endpoint","present")))

offset <- 0.12

emm_time <- emmeans(m_time, ~ Group) %>% as.data.frame() %>%
  mutate(Group = factor(Group, levels = c("endpoint","present")))

dat_pos_t <- dat %>%
  mutate(
    x_box = as.numeric(Group) + ifelse(Group == "present", +offset, -offset),
    x_emm = as.numeric(Group) + ifelse(Group == "present", -offset, +offset)
  )
emm_pos_t <- emm_time %>%
  mutate(x_emm = as.numeric(Group) + ifelse(Group == "present", -offset, +offset))

outliers_t <- dat_pos_t %>%
  group_by(Group) %>%
  mutate(is_out = flag_outliers(TenPercp)) %>%
  ungroup() %>%
  filter(is_out)

p_time <- ggplot() +
  geom_boxplot(
    data = dat_pos_t,
    aes(x = x_box, y = TenPercp, colour = Group, fill = Group, group = Group),
    width = 0.30, linewidth = 0.35, alpha = 0.15, outlier.shape = NA
  ) +
  geom_point(
    data = outliers_t,
    aes(x = x_box, y = TenPercp),
    shape = 21, size = 0.5, stroke = 0.35, fill = "black", colour = "black"
  ) +
  geom_errorbar(
    data = emm_pos_t,
    aes(x = x_emm, y = emmean, ymin = lower.CL, ymax = upper.CL, colour = Group),
    width = 0.06, linewidth = 0.45
  ) +
  geom_point(
    data = emm_pos_t,
    aes(x = x_emm, y = emmean, colour = Group, shape = Group),
    fill = "white", size = 1.0, stroke = 0.7
  ) +
  geom_line(
    data = emm_pos_t,
    aes(x = x_emm, y = emmean, group = 1),
    colour = "black", linetype = "dashed", linewidth = 0.4
  ) +
  scale_x_continuous(breaks = c(1, 2), labels = c("endpoint","present"), expand = c(0.2, 0.2)) +
  scale_shape_manual(values = c("endpoint" = 21, "present" = 24)) +
  scale_colour_manual(values = group_colors) +
  scale_fill_manual(values   = group_colors) +
  scale_y_continuous(limits = c(1, 9), breaks = seq(1, 9, 2), expand = c(0, 0)) +
  labs(x = NULL, y = "Rating") +
  theme(
    axis.line   = element_line(colour = "black", linewidth = 0.35),
    axis.title  = element_text(size = 7, colour = "black"),
    axis.text   = element_text(size = 7, colour = "black"),
    legend.position = "none",
    plot.margin = ggplot2::margin(3, 3, 2, 2)
  ) +
  patchwork::plot_annotation(
    title = "(b) Time perception",
    theme = theme(plot.title = element_text(size = 7, colour = "black", face = "bold",
                                            margin = ggplot2::margin(b = 2)))
  )

ggsave("pic_time.pdf", plot = p_time, width = 2, height = 2, units = "in", device = cairo_pdf)

####################################################################################################
### 5) Resource allocation strategy ----------------------------------------------------------------
####################################################################################################

# Difference-score ANCOVAs
dfc <- dat %>%
  transmute(
    SubjNum,
    Group  = relevel(factor(Group), ref = "present"),
    Gender = factor(Gender, levels = c("male","female")),
    Age_z,
    spend  = Account_Immediate,
    short  = Account_ShortTerm,
    long   = Account_LongTerm
  ) %>%
  mutate(
    D_S  = short - spend,   # short vs spending
    D_L  = long  - spend,   # long  vs spending
    D_LS = long  - short    # long  vs short
  )

# D_S
m_DS     <- lm(D_S  ~ Group + Age_z + Gender, data = dfc)
a3_DS    <- withr::with_options(list(contrasts = c("contr.sum","contr.poly")),
                                car::Anova(m_DS, type = 3))
eta_DS   <- eta_from_Anova(as.data.frame(a3_DS), "Group")
emm_DS   <- emmeans(m_DS, ~ Group)
cmp_DS   <- pairs(emm_DS, infer = c(TRUE, TRUE))

# D_L
m_DL     <- lm(D_L  ~ Group + Age_z + Gender, data = dfc)
a3_DL    <- withr::with_options(list(contrasts = c("contr.sum","contr.poly")),
                                car::Anova(m_DL, type = 3))
eta_DL   <- eta_from_Anova(as.data.frame(a3_DL), "Group")
emm_DL   <- emmeans(m_DL, ~ Group)
cmp_DL   <- pairs(emm_DL, infer = c(TRUE, TRUE))

# D_LS
m_DLS    <- lm(D_LS ~ Group + Age_z + Gender, data = dfc)
a3_DLS   <- withr::with_options(list(contrasts = c("contr.sum","contr.poly")),
                                car::Anova(m_DLS, type = 3))
eta_DLS  <- eta_from_Anova(as.data.frame(a3_DLS), "Group")
emm_DLS  <- emmeans(m_DLS, ~ Group)
cmp_DLS  <- pairs(emm_DLS, infer = c(TRUE, TRUE))

# Effect size calculation (saving vs. spending)
b0_DS <- coef(summary(m_DS))["(Intercept)", "Estimate"]
b0_DL <- coef(summary(m_DL))["(Intercept)", "Estimate"]

a3_DS  <- withr::with_options(list(contrasts=c("contr.sum","contr.poly")), car::Anova(m_DS, type=3))
a3_DL  <- withr::with_options(list(contrasts=c("contr.sum","contr.poly")), car::Anova(m_DL, type=3))

F_DS_i <- a3_DS["(Intercept)", "F value"];  df2_DS <- a3_DS["Residuals", "Df"]
F_DL_i <- a3_DL["(Intercept)", "F value"];  df2_DL <- a3_DL["Residuals", "Df"]

eta_DS_i <- F_DS_i / (F_DS_i + df2_DS)
eta_DL_i <- F_DL_i / (F_DL_i + df2_DL)

# Multiplicity control across 3 group contrasts (Holm)
tab_cmp <- bind_rows(
  cbind(measure = "D_S (short - spending)", as.data.frame(cmp_DS)),
  cbind(measure = "D_L (long  - spending)", as.data.frame(cmp_DL)),
  cbind(measure = "D_LS (long - short)",    as.data.frame(cmp_DLS))
) %>%
  mutate(padj_holm = p.adjust(p.value, method = "holm"))

# Compact ANCOVA summary (F, df, p, partial ¦Ç2)
summ_line <- function(a3, eta_obj, label) {
  data.frame(
    measure = label,
    F       = as.numeric(a3["Group", "F value"]),
    df1     = as.numeric(a3["Group", "Df"]),
    df2     = as.numeric(a3["Residuals", "Df"]),
    p       = as.numeric(a3["Group", "Pr(>F)"]),
    eta_p2  = as.numeric(eta_obj$Eta2_partial),
    CI_low  = if (!is.null(eta_obj$CI_low))  as.numeric(eta_obj$CI_low)  else NA_real_,
    CI_high = if (!is.null(eta_obj$CI_high)) as.numeric(eta_obj$CI_high) else NA_real_
  )
}

tab_ancova <- bind_rows(
  summ_line(as.data.frame(a3_DS),  eta_DS,  "D_S (short - spending)"),
  summ_line(as.data.frame(a3_DL),  eta_DL,  "D_L (long  - spending)"),
  summ_line(as.data.frame(a3_DLS), eta_DLS, "D_LS (long - short)")
)

print(a3_DS); print(a3_DL); print(a3_DLS)  # Type-III tables
tab_ancova                                # For reporting
tab_cmp                                    # Holm-adjusted pairwise (Group)

# ROBUSTNESS: short vs long only (same sample as D_LS)
dat_SL <- dat %>%
  transmute(
    SubjNum,
    Group  = relevel(factor(Group), ref = "present"),
    Gender = factor(Gender, levels = c("male","female")),
    Age_z,
    short  = Account_ShortTerm,
    long   = Account_LongTerm
  ) %>%
  filter(!is.na(short) & !is.na(long)) %>%
  tidyr::pivot_longer(c(short, long), names_to = "Account", values_to = "prop") %>%
  mutate(Account = factor(Account, levels = c("short","long")))

m_SL <- lmer(prop ~ Group*Account + Age_z + Gender + (1|SubjNum), data = dat_SL)
withr::with_options(
  list(contrasts = c("contr.sum","contr.poly")),
  print(car::Anova(m_SL, type = 3))  # Group:Account ¡Ö p .06
)

emm_SL  <- emmeans(m_SL, ~ Account | Group)
print(summary(contrast(emm_SL, list("Long - Short" = c(-1, 1)), by = "Group"),
              infer = c(TRUE, TRUE)))

emm_SL2 <- emmeans(m_SL, ~ Account * Group)
print(summary(contrast(emm_SL2,
                       list("Delta(Long-Short): endpoint - present" = c(-1, 1, 1, -1))),
              infer = c(TRUE, TRUE)))

# Plot
# Adjusted means (¡À95% CI) for each Account ¡Á Group
em_SL_plot <- emmeans(m_SL, ~ Account * Group) %>%
  as.data.frame() %>%
  mutate(
    Group   = factor(Group,  levels = c("endpoint","present")),
    Account = factor(Account, levels = c("short","long"))
  )

# Aesthetics setting
group_colors <- c("present" = "#4169E1", "endpoint" = "#B22222")
group_shapes <- c("endpoint" = 21, "present" = 24)  # endpoint=¡ð, present=¡÷
pd <- position_dodge(width = 0.25)

p_SL <- ggplot(em_SL_plot,
               aes(x = Account, y = emmean,
                   group = Group, colour = Group, shape = Group)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.10, size = 0.45, position = pd) +
  geom_point(position = pd, fill = "white", size = 1, stroke = 0.7) +
  geom_line(position = pd, linewidth = 0.45) +
  scale_colour_manual(values = group_colors) +
  scale_shape_manual(values  = group_shapes) +
  scale_y_continuous(limits = c(0, 1),  oob = scales::squish,
                     breaks = seq(0, 1, by = 0.25),
                     expand = c(0, 0),
                     labels = scales::label_percent(scale = 100)) +
  labs(x = NULL, y = "Allocation proportion (%)") +
  theme_classic(base_size = 8) +
  theme(
    axis.line        = element_line(colour = "black", linewidth = 0.35),
    axis.title       = element_text(size = 7, colour = "black"),
    axis.text        = element_text(size = 7, colour = "black"),
    panel.background = element_rect(fill = "transparent"),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    legend.position  = "none",
    plot.margin      = ggplot2::margin(3, 3, 2, 2)
  ) +
  plot_annotation(
    title = "(c) Resource allocation",
    theme = theme(plot.title = element_text(size = 7, colour = "black",
                                            face = "bold",
                                            margin = ggplot2::margin(b = 2)))
  )

ggsave("pic_ra.pdf", plot = p_SL, width = 2, height = 2,
       units = "in", device = cairo_pdf)

####################################################################################################
# END
####################################################################################################
