# Meditation: Endpoint-focused vs. Present-focused
# Behavioral data analysis & plotting
# For: mental image
# Requires: "rawdata.xlsx" (sheet = "rawdata2")
# Programmed by Feng XIAO (updated on 2026-3-19)

####################################################################################################
### 0) Preparation ---------------------------------------------------------------------------------
####################################################################################################

# Load packages actually used
pkg <- c(
  "readxl", "dplyr", "car", "emmeans", "effectsize"
)
lapply(pkg, require, character.only = TRUE)

# Working directory (RStudio only)
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Function for within-group analysis summary
get_baseline_stats <- function(model, baseline_value) {
  s  <- summary(model)$coefficients
  lh <- car::linearHypothesis(model, paste0("(Intercept) = ", baseline_value))
  
  data.frame(
    estimate = s["(Intercept)", "Estimate"],
    se       = s["(Intercept)", "Std. Error"],
    F        = lh[2, "F"],
    df1      = lh[2, "Df"],
    df2      = lh[2, "Res.Df"],   
    p        = lh[2, "Pr(>F)"]
  )
}

# Function for partial eta squared
get_partial_eta2 <- function(model, hypothesis) {
  lh <- car::linearHypothesis(model, hypothesis)
  
  F_value <- lh[2, "F"]
  df1     <- lh[2, "Df"]
  df2     <- lh[2, "Res.Df"]
  
  eta_p2 <- (F_value * df1) / (F_value * df1 + df2)
  return(eta_p2)
}

####################################################################################################
### 1) Load data -----------------------------------------------------------------------------------
####################################################################################################

dat <- read_excel("rawdata.xlsx", sheet = "rawdata2")

dat <- dat %>%
  mutate(
    Age_z = scale(Age)[,1],
    Gender = factor(Gender),
    Group = factor(Group) 
  )

dat_present  <- dat %>% filter(Group == "present")
dat_endpoint <- dat %>% filter(Group == "endpoint")

####################################################################################################
### 2) Present-focused ---------------------------------------------------------------------------
####################################################################################################

m_ev_p   <- lm(EmotionalValence ~ Age_z + Gender, data = dat_present)
m_self_p <- lm(SelfRelatedness  ~ Age_z + Gender, data = dat_present)
m_soc_p  <- lm(SocialCloseness  ~ Age_z + Gender, data = dat_present)

# p
p_present <- c(
  EmotionalValence = car::linearHypothesis(m_ev_p,   "(Intercept) = 5")[2, "Pr(>F)"],
  SelfRelatedness  = car::linearHypothesis(m_self_p, "(Intercept) = 0")[2, "Pr(>F)"],
  SocialCloseness  = car::linearHypothesis(m_soc_p,  "(Intercept) = 0")[2, "Pr(>F)"]
)

# partial eta squared
eta2p_present <- c(
  EmotionalValence = get_partial_eta2(m_ev_p,   "(Intercept) = 5"),
  SelfRelatedness  = get_partial_eta2(m_self_p, "(Intercept) = 0"),
  SocialCloseness  = get_partial_eta2(m_soc_p,  "(Intercept) = 0")
)

stats_present <- rbind(
  EmotionalValence = get_baseline_stats(m_ev_p,   5),
  SelfRelatedness  = get_baseline_stats(m_self_p, 0),
  SocialCloseness  = get_baseline_stats(m_soc_p,  0)
)

stats_present$p_fdr <- p.adjust(stats_present$p, method = "fdr")
stats_present$eta_p2 <- eta2p_present

cat("\n--- Present group baseline-aware stats ---\n")
print(stats_present)

####################################################################################################
### 3) Endpoint-focused ----------------------------------------------------------------------------------------
####################################################################################################

m_ev_e   <- lm(EmotionalValence ~ Age_z + Gender, data = dat_endpoint)
m_self_e <- lm(SelfRelatedness  ~ Age_z + Gender, data = dat_endpoint)
m_soc_e  <- lm(SocialCloseness  ~ Age_z + Gender, data = dat_endpoint)

# p
p_endpoint <- c(
  EmotionalValence = car::linearHypothesis(m_ev_e,   "(Intercept) = 5")[2, "Pr(>F)"],
  SelfRelatedness  = car::linearHypothesis(m_self_e, "(Intercept) = 0")[2, "Pr(>F)"],
  SocialCloseness  = car::linearHypothesis(m_soc_e,  "(Intercept) = 0")[2, "Pr(>F)"]
)

# partial eta squared
eta2p_endpoint <- c(
  EmotionalValence = get_partial_eta2(m_ev_e,   "(Intercept) = 5"),
  SelfRelatedness  = get_partial_eta2(m_self_e, "(Intercept) = 0"),
  SocialCloseness  = get_partial_eta2(m_soc_e,  "(Intercept) = 0")
)

stats_endpoint <- rbind(
  EmotionalValence = get_baseline_stats(m_ev_e,   5),
  SelfRelatedness  = get_baseline_stats(m_self_e, 0),
  SocialCloseness  = get_baseline_stats(m_soc_e,  0)
)

stats_endpoint$p_fdr <- p.adjust(stats_endpoint$p, method = "fdr")
stats_endpoint$eta_p2 <- eta2p_endpoint

cat("\n--- Endpoint group baseline-aware stats ---\n")
print(stats_endpoint)

####################################################################################################
### 4) Between-group (Endpoint vs Present) ----------------------------------------------------------
####################################################################################################

m_ev_g   <- lm(EmotionalValence ~ Group + Age_z + Gender, data = dat)
m_self_g <- lm(SelfRelatedness  ~ Group + Age_z + Gender, data = dat)
m_soc_g  <- lm(SocialCloseness  ~ Group + Age_z + Gender, data = dat)

# Group main effect p
p_group <- c(
  EmotionalValence = car::Anova(m_ev_g,   type = 3)["Group", "Pr(>F)"],
  SelfRelatedness  = car::Anova(m_self_g, type = 3)["Group", "Pr(>F)"],
  SocialCloseness  = car::Anova(m_soc_g,  type = 3)["Group", "Pr(>F)"]
)

# FDR
p_group_fdr <- p.adjust(p_group, method = "fdr")

# F values
F_group <- c(
  EmotionalValence = car::Anova(m_ev_g,   type = 3)["Group", "F value"],
  SelfRelatedness  = car::Anova(m_self_g, type = 3)["Group", "F value"],
  SocialCloseness  = car::Anova(m_soc_g,  type = 3)["Group", "F value"]
)

# df
df_group <- c(
  EmotionalValence = car::Anova(m_ev_g,   type = 3)["Group", "Df"],
  SelfRelatedness  = car::Anova(m_self_g, type = 3)["Group", "Df"],
  SocialCloseness  = car::Anova(m_soc_g,  type = 3)["Group", "Df"]
)

# partial eta^2
eta_group <- c(
  EmotionalValence = effectsize::eta_squared(m_ev_g, partial = TRUE)$Eta2_partial[1],
  SelfRelatedness  = effectsize::eta_squared(m_self_g, partial = TRUE)$Eta2_partial[1],
  SocialCloseness  = effectsize::eta_squared(m_soc_g, partial = TRUE)$Eta2_partial[1]
)

stats_group <- data.frame(
  F       = F_group,
  df      = df_group,
  p       = p_group,
  p_fdr   = p_group_fdr,
  eta_p2  = eta_group
)

cat("\n--- Between-group effects (Group) ---\n")
print(stats_group)

####################################################################################################
# END
####################################################################################################