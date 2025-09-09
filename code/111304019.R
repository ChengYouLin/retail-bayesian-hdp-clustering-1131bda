# 1131_貝氏資料分析_Final
# 111304019 統計三 林承佑

# 載入套件
library(rjags)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(showtext)
showtext_auto()



# 讀取資料和處理資料格式

data <- read.csv("111304019.csv")
data <- data[, -c(1)]
data <- data %>%
  mutate(
    Quarter = as.numeric(factor(季度, levels = c("第一季", "第二季", "第三季", "第四季"))),
    County = as.numeric(as.factor(縣市)),
    RegionType = as.numeric(factor(地區, levels = c("一般縣市", "直轄市"))),
    ProductType = as.numeric(factor(品項定位, levels = c("中性", "特定性別"))),
    Style = as.numeric(factor(款式, levels = c("一般款", "自動扣", "針棒")))
  )

# 確認每個因子變數的 Level 與對應編號
quarter_levels <- levels(factor(data$季度, levels = c("第一季", "第二季", "第三季", "第四季")))
quarter_mapping <- setNames(seq_along(quarter_levels), quarter_levels)
county_levels <- levels(as.factor(data$縣市))
county_mapping <- setNames(seq_along(county_levels), county_levels)
region_levels <- levels(factor(data$地區, levels = c("一般縣市", "直轄市")))
region_mapping <- setNames(seq_along(region_levels), region_levels)
product_levels <- levels(factor(data$品項定位, levels = c("中性", "特定性別")))
product_mapping <- setNames(seq_along(product_levels), product_levels)
style_levels <- levels(factor(data$款式, levels = c("一般款", "自動扣", "針棒")))
style_mapping <- setNames(seq_along(style_levels), style_levels)

list(
  Quarter = quarter_mapping,
  County = county_mapping,
  RegionType = region_mapping,
  ProductType = product_mapping,
  Style = style_mapping
)


# 簡單查看銷售分佈
sales_by_style <- data %>%
  group_by(Style) %>%
  summarize(MeanSales = mean(銷售總金額), SD = sd(銷售總金額))

ggplot(data, aes(x = Style, y = 銷售總金額)) +
  geom_boxplot() +
  ggtitle("Sales Distribution by Style") +
  theme_minimal()



#########模型準備#########

# 準備給jags 的資料

jags_data <- list(
  N = nrow(data),
  Sales = data$Sales,
  Quarter = data$Quarter,
  RegionType = data$RegionType,
  ProductType = data$ProductType,
  Style = data$Style
)

# 初始for MCMC
inits <- list(
  list(beta_0 = 0, sigma = 1, df = 5),
  list(beta_0 = 1, sigma = 2, df = 10),
  list(beta_0 = -1, sigma = 1.5, df = 3)
)




################
## 第一個模型 ##
################

jags_model <- "
model {
  for (i in 1:N) {
    Sales[i] ~ dt(mu[i], tau, df)  # t-distribution for robustness
    mu[i] <- beta_0 + beta_quarter[Quarter[i]] + beta_region[RegionType[i]] + \
             beta_product[ProductType[i]] + beta_style[Style[i]] + \
             beta_interaction[ProductType[i], Style[i]]
  }

  # Priors for fixed effects
  beta_0 ~ dnorm(0, 0.001)

  for (j in 1:4) {
    beta_quarter[j] ~ dnorm(0, 0.001)
  }

  for (k in 1:2) {
    beta_region[k] ~ dnorm(0, 0.001)
  }

  for (m in 1:2) {
    beta_product[m] ~ dnorm(0, 0.001)
  }

  for (n in 1:3) {
    beta_style[n] ~ dnorm(0, 0.001)
  }

  for (p in 1:2) {
    for (q in 1:3) {
      beta_interaction[p, q] ~ dnorm(0, 0.001)
    }
  }

  # Priors for variance components
  tau <- 1 / (sigma^2)
  sigma ~ dunif(0, 100)
  df ~ dexp(1)  # Degrees of freedom for t-distribution
}
"

# 執行
model <- jags.model(textConnection(jags_model), data = jags_data, inits = inits, n.chains = 3, n.adapt = 10000)
update(model, 5000)
samples <- coda.samples(model, variable.names = c("beta_0", "beta_quarter", "beta_region", "beta_product", "beta_style", "beta_interaction"), n.iter = 5000)
summary(samples)


mcmc_samples <- do.call(rbind, lapply(samples, as.data.frame))
mcmc_samples <- as.data.frame(mcmc_samples)

# 各項參數貢獻值
effect_summary <- mcmc_samples %>%
  summarize(across(starts_with("beta"), mean)) %>%
  pivot_longer(cols = everything(), names_to = "Factor", values_to = "Effect")

ggplot(effect_summary, aes(x = Factor, y = abs(Effect))) +
  geom_bar(stat = "identity") +
  ggtitle("Effect Size by Factor") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##############
# likelihood #
##############

## 款式
# 款式的均值和精度參數
style_effects  <- mcmc_samples %>%
  dplyr::select(starts_with("beta_style")) %>%
  pivot_longer(cols = everything(), names_to = "Style", values_to = "Effect")

# 計算每款式的最大似然估計值
style_mle <- style_effects %>%
  group_by(Style) %>%
  summarize(MLE = mean(Effect), CI_Lower = quantile(Effect, 0.025), CI_Upper = quantile(Effect, 0.975))
print(style_mle)

ggplot(style_effects, aes(x = Effect, fill = Style)) +
  geom_density(alpha = 0.5) +
  ggtitle("Posterior (Likelihood) Distributions for Each Style") +
  theme_minimal()


## 品項定位
# 品項定位的均值和精度參數
style_effects <- mcmc_samples %>%
  dplyr::select(starts_with("beta_product")) %>%
  pivot_longer(cols = everything(), names_to = "product", values_to = "Effect")

# 計算每款式的最大似然估計值
style_pro <- style_effects %>%
  group_by(product) %>%
  summarize(MLE = mean(Effect), CI_Lower = quantile(Effect, 0.025), CI_Upper = quantile(Effect, 0.975))
print(style_pro)

ggplot(style_effects, aes(x = Effect, fill = product)) +
  geom_density(alpha = 0.5) +
  ggtitle("Posterior (Likelihood) Distributions for Each Style") +
  theme_minimal()

## 地區
# 提取地區效應的後驗分布
region_effects <- mcmc_samples %>%
  dplyr::select(starts_with("beta_region")) %>%
  pivot_longer(cols = everything(), names_to = "Region", values_to = "Effect")

# 計算每地區的最大似然估計值
region_mle <- region_effects %>%
  group_by(Region) %>%
  summarize(MLE = mean(Effect), CI_Lower = quantile(Effect, 0.025), CI_Upper = quantile(Effect, 0.975))
print(region_mle)

# 後驗分布作為似然的近似
ggplot(region_effects, aes(x = Effect, fill = Region)) +
  geom_density(alpha = 0.5) +
  ggtitle("Posterior (Likelihood) Distributions for Each Region") +
  theme_minimal()

# 觀察季度這項參數，他對於銷售總金額的狀況
quarter_effects <- mcmc_samples %>%
  dplyr::select(starts_with("beta_quarter")) %>%
  pivot_longer(cols = everything(), names_to = "Quarter", values_to = "Effect")

ggplot(quarter_effects, aes(x = Quarter, y = Effect)) +
  geom_boxplot() +
  ggtitle("Posterior Distribution of Quarterly Effects") +
  theme_minimal()

# 觀察不同組合下的likelihood
interaction_effects <- mcmc_samples %>%
  dplyr::select(starts_with("beta_interaction")) %>%
  pivot_longer(cols = everything(), names_to = "Interaction", values_to = "Effect")

ggplot(interaction_effects, aes(x = Interaction, y = Effect)) +
  geom_boxplot() +
  ggtitle("Posterior Distribution of Interaction Effects") +
  theme_minimal()


# 觀察 款式 這項參數，他對於銷售總金額的狀況
style_effects <- mcmc_samples %>%
  dplyr::select(starts_with("beta_style")) %>%
  pivot_longer(cols = everything(), names_to = "Style", values_to = "Effect")

ggplot(style_effects, aes(x = Style, y = Effect)) +
  geom_boxplot() +
  ggtitle("Posterior Distribution of Style Effects") +
  theme_minimal()


#比較季度間的和不同產品組合的差異
interaction_quarterly_effects <- mcmc_samples %>%
  dplyr::select(starts_with("beta_interaction")) %>%
  pivot_longer(cols = everything(), names_to = "Combination", values_to = "Effect") %>%
  mutate(Quarter = factor(rep(1:4, length.out = n()), labels = c("Q1", "Q2", "Q3", "Q4")))

# 繪製密度分布圖，分季度 -> 每季度一張圖
ggplot(interaction_quarterly_effects, aes(x = Effect, fill = Combination)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Quarter, ncol = 2) +  
  ggtitle("Density Plots of Interaction Effects by Quarter") +
  xlab("Interaction Effect") +
  ylab("Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 熵值計算
entropy_quarterly <- interaction_quarterly_effects %>%
  group_by(Quarter, Combination) %>%
  summarize(Entropy = -sum(log(Effect^2)))

# 使用點線圖展示熵值隨季度的變化
ggplot(entropy_quarterly, aes(x = Quarter, y = Entropy, group = Combination, color = Combination)) +
  geom_line(size = 1) +             
  geom_point(size = 3) +               
  scale_x_discrete(labels = c("Q1", "Q2", "Q3", "Q4")) + 
  ggtitle("Entropy Trends Across Quarters for Each Combination") +
  xlab("Quarter") +
  ylab("Entropy") +
  theme_minimal() +
  theme(legend.position = "bottom")  

################
## 第二個模型 ##
################

#########模型準備#########

# 準備給jags 的資料
data$Sales <- data$銷售總金額

jags_data <- list(
  N = nrow(data),
  Sales = data$Sales,
  Quarter = data$Quarter,
  RegionType = data$RegionType,
  ProductType = data$ProductType,
  Style = data$Style,
  mean_quarter = tapply(data$Sales, data$Quarter, mean, na.rm = TRUE),
  prec_quarter = 1 / tapply(data$Sales, data$Quarter, var, na.rm = TRUE),
  mean_region = tapply(data$Sales, data$RegionType, mean, na.rm = TRUE),
  prec_region = 1 / tapply(data$Sales, data$RegionType, var, na.rm = TRUE),
  mean_product = tapply(data$Sales, data$ProductType, mean, na.rm = TRUE),
  prec_product = 1 / tapply(data$Sales, data$ProductType, var, na.rm = TRUE),
  mean_style = tapply(data$Sales, data$Style, mean, na.rm = TRUE),
  prec_style = 1 / tapply(data$Sales, data$Style, var, na.rm = TRUE)
)



jags_model_quarter_adjusted <- "
model {
  for (i in 1:N) {
    Sales[i] ~ dt(mu[i], tau, df)  # t-distribution for robustness
    mu[i] <- beta_0 + beta_quarter + beta_region[RegionType[i]] + 
             beta_product[ProductType[i]] + beta_style[Style[i]] + 
             beta_interaction[ProductType[i], Style[i]]
  }

  # Priors for fixed effects
  beta_0 ~ dnorm(0, 0.001)
  beta_quarter ~ dnorm(0, 0.001)  # Single prior for all quarters

  for (k in 1:2) {
    beta_region[k] ~ dnorm(0, 0.001)
  }

  for (m in 1:2) {
    beta_product[m] ~ dnorm(0, 0.001)
  }

  for (n in 1:3) {
    beta_style[n] ~ dnorm(0, 0.001)
  }

  for (p in 1:2) {
    for (q in 1:3) {
      beta_interaction[p, q] ~ dnorm(0, 0.001)
    }
  }

  # Priors for variance components
  tau <- 1 / (sigma^2)
  sigma ~ dunif(0, 100)
  df ~ dexp(1)  # Degrees of freedom for t-distribution
}
"

model <- jags.model(textConnection(jags_model_quarter_adjusted), data = jags_data, n.chains = 3, n.adapt = 1000)
update(model, 1000) 
samples <- coda.samples(model, variable.names = c("beta_0", "beta_quarter", "beta_region", "beta_interaction", "z", "beta_style","beta_product"), n.iter = 5000)

mcmc_samples <- do.call(rbind, lapply(samples, as.data.frame))
mcmc_samples <- as.data.frame(mcmc_samples)

# 後驗分布
beta_quarter_posterior <- mcmc_samples %>% 
  dplyr::select(starts_with("beta_quarter")) %>%
  rename(Effect = beta_quarter)
summary(beta_quarter_posterior)

# 計算均值與置信區間
mean_effect <- mean(beta_quarter_posterior$Effect)
ci_lower <- quantile(beta_quarter_posterior$Effect, 0.025)
ci_upper <- quantile(beta_quarter_posterior$Effect, 0.975)
cat("Mean:", mean_effect, "95% CI:", ci_lower, "-", ci_upper, "\n")

ggplot(beta_quarter_posterior, aes(x = Effect)) +
  geom_density(fill = "blue", alpha = 0.5) +
  ggtitle("Posterior Distribution of Beta_Quarter") +
  xlab("Quarter Effect") +
  ylab("Density") +
  theme_minimal()

ggplot(beta_quarter_posterior, aes(x = "Quarter Effect", y = Effect)) +
  geom_boxplot() +
  ggtitle("Boxplot of Posterior Beta_Quarter") +
  xlab("") +
  ylab("Effect") +
  theme_minimal()


## 款式
# 款式的均值和精度參數
style_effects <- mcmc_samples %>%
  dplyr::select(starts_with("beta_style")) %>%
  pivot_longer(cols = everything(), names_to = "Style", values_to = "Effect")

# 計算每款式的最大似然估計值
print("原本的：")
print(style_mle)
style_mle <- style_effects %>%
  group_by(Style) %>%
  summarize(MLE = mean(Effect), CI_Lower = quantile(Effect, 0.025), CI_Upper = quantile(Effect, 0.975))
print(style_mle)

ggplot(style_effects, aes(x = Effect, fill = Style)) +
  geom_density(alpha = 0.5) +
  ggtitle("Posterior (Likelihood) Distributions for Each Style") +
  theme_minimal()


## 品項定位
# 品項定位的均值和精度參數
style_effects <- mcmc_samples %>%
  dplyr::select(starts_with("beta_product")) %>%
  pivot_longer(cols = everything(), names_to = "product", values_to = "Effect")

# 計算每款式的最大似然估計值
print("原本的：")
print(style_pro)
style_pro <- style_effects %>%
  group_by(product) %>%
  summarize(MLE = mean(Effect), CI_Lower = quantile(Effect, 0.025), CI_Upper = quantile(Effect, 0.975))
print(style_pro)

ggplot(style_effects, aes(x = Effect, fill = product)) +
  geom_density(alpha = 0.5) +
  ggtitle("Posterior (Likelihood) Distributions for Each Style") +
  theme_minimal()



################
## 第三個模型 ##
################

jags_model_region_prior <- "
model {
  for (i in 1:N) {
    Sales[i] ~ dt(mu[i], tau, df)
    mu[i] <- beta_0 + beta_quarter[Quarter[i]] + beta_region[RegionType[i]] +
             beta_product[ProductType[i]] + beta_style[Style[i]] +
             beta_interaction[ProductType[i], Style[i]]
  }

  # Priors for fixed effects
  beta_0 ~ dnorm(0, 0.001)

  for (j in 1:4) {
    beta_quarter[j] ~ dnorm(0, 0.001)
  }

  # Differentiated priors for regions
  beta_region[1] ~ dnorm(0, 0.002)  # Higher variance for general regions
  beta_region[2] ~ dnorm(0, 0.0005)  # Lower variance for metropolitan areas

  for (m in 1:2) {
    beta_product[m] ~ dnorm(0, 0.001)
  }

  for (n in 1:3) {
    beta_style[n] ~ dnorm(0, 0.001)
  }

  for (p in 1:2) {
    for (q in 1:3) {
      beta_interaction[p, q] ~ dnorm(0, 0.001)
    }
  }

  # Priors for variance components
  tau <- 1 / (sigma^2)
  sigma ~ dunif(0, 100)
  df ~ dexp(1)  # Degrees of freedom for t-distribution
}
"

model <- jags.model(textConnection(jags_model_region_prior), data = jags_data, n.chains = 3, n.adapt = 1000)
update(model, 1000) 
samples <- coda.samples(model, variable.names = c("beta_0", "beta_quarter", "beta_region", "beta_interaction", "z"), n.iter = 5000)


mcmc_samples <- do.call(rbind, lapply(samples, as.data.frame))
mcmc_samples <- as.data.frame(mcmc_samples)
posterior_samples <- mcmc_samples

ggplot(region_effects, aes(x = Region, y = Effect, fill = Region)) +
  geom_boxplot() +
  ggtitle("Boxplot of Region Effects") +
  ylab("Effect") +
  theme_minimal()

region_effects <- posterior_samples %>%
  dplyr::select(starts_with("beta_region")) %>%
  pivot_longer(cols = everything(), names_to = "Region", values_to = "Effect") %>%
  mutate(
    Region = case_when(
      str_detect(Region, "beta_region.1") ~ "一般縣市",
      str_detect(Region, "beta_region.2") ~ "直轄市"
    )
  )

ggplot(region_effects, aes(x = Effect, fill = Region)) +
  geom_density(alpha = 0.5) +
  ggtitle("Posterior Distributions of Region Effects") +
  xlab("Effect") +
  ylab("Density") +
  theme_minimal()


region_summary <- region_effects %>%
  group_by(Region) %>%
  summarize(
    Mean_Effect = mean(Effect),
    Lower_CI = quantile(Effect, 0.025),
    Upper_CI = quantile(Effect, 0.975),
    .groups = "drop"
  )
print(region_summary)


################
## 第四個模型 ##
################
jags_model_latent_mixture <- "
model {
  for (i in 1:N) {
    Sales[i] ~ dt(mu[i], tau, df)
    mu[i] <- beta_0 + beta_quarter[Quarter[i]] + beta_region[RegionType[i]] +
             beta_product[ProductType[i]] + beta_style[Style[i]] +
             z[i] * beta_interaction[ProductType[i], Style[i]]  # Latent group effect
  }

  # Priors for fixed effects
  beta_0 ~ dnorm(0, 0.001)

  for (j in 1:4) {
    beta_quarter[j] ~ dnorm(0, 0.001)
  }

  for (k in 1:2) {
    beta_region[k] ~ dnorm(0, 0.001)
  }

  for (m in 1:2) {
    beta_product[m] ~ dnorm(0, 0.001)
  }

  for (n in 1:3) {
    beta_style[n] ~ dnorm(0, 0.001)
  }

  for (p in 1:2) {
    for (q in 1:3) {
      beta_interaction[p, q] ~ dnorm(0, 0.001)
    }
  }

  # Latent mixture model
  for (i in 1:N) {
    z[i] ~ dbern(pi[i])  # Binary latent variable
    logit(pi[i]) <- alpha + beta_latent[RegionType[i]]
  }

  alpha ~ dnorm(0, 0.001)
  beta_latent[1] ~ dnorm(0, 0.001)
  beta_latent[2] ~ dnorm(0, 0.001)

  # Priors for variance components
  tau <- 1 / (sigma^2)
  sigma ~ dunif(0, 100)
  df ~ dexp(1)
}
"

model <- jags.model(textConnection(jags_model_latent_mixture), data = jags_data, n.chains = 3, n.adapt = 1000)
update(model, 1000) 
samples <- coda.samples(model, variable.names = c("beta_0", "beta_quarter", "beta_region", "beta_interaction", "z"), n.iter = 5000)


samples_df <- do.call(rbind, lapply(samples, as.data.frame))
samples_df <- as.data.frame(samples_df)
posterior_samples <- mcmc_samples


# z[i] 的後驗均值（分組概率）
latent_probs <- samples_df %>%
  select(starts_with("z[")) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything(), names_to = "DataPoint", values_to = "PosteriorMean")

# 分組判斷：將概率大於 ？？ 的分配到組 1，否則分配到組 0
# 需要手動修改分組閾值
latent_probs <- latent_probs %>%
  mutate(Group = ifelse(PosteriorMean > 0.9, 1, 0))

group_summary <- latent_probs %>%
  group_by(Group) %>%
  summarise(Count = n())
print(group_summary)

# 合併潛在分組
data_with_group <- cbind(data, Group = latent_probs$Group)

# 按群統計特徵
group_characteristics <- data_with_group %>%
  group_by(Group, RegionType, ProductType, Style) %>%
  summarise(AverageSales = mean(Sales), Count = n(), .groups = "drop")
print(group_characteristics)

ggplot(latent_probs, aes(x = PosteriorMean, fill = as.factor(Group))) +
  geom_density(alpha = 0.5) +
  ggtitle("Posterior Distribution of Latent Variables (z[i])") +
  theme_minimal() +
  xlab("Posterior Mean of z[i]") +
  ylab("Density") +
  scale_fill_discrete(name = "Group")

# 每組的後驗概率特徵
group_posterior_stats <- latent_probs %>%
  group_by(Group) %>%
  summarise(AveragePosterior = mean(PosteriorMean),
            SDPosterior = sd(PosteriorMean),
            MinPosterior = min(PosteriorMean),
            MaxPosterior = max(PosteriorMean))
print(group_posterior_stats)

# 比較組間銷售額的顯著性
sales_by_group <- data_with_group %>%
  group_by(Group) %>%
  summarise(AverageSales = mean(Sales), SD = sd(Sales), Count = n())

# t 檢驗檢查組間差異
t_test_result <- t.test(Sales ~ Group, data = data_with_group)
print(t_test_result)


# 分組內交互特徵分布
interaction_analysis <- data_with_group %>%
  group_by(Group, RegionType, ProductType, Style) %>%
  summarise(AverageSales = mean(Sales), Count = n(), .groups = "drop")

ggplot(interaction_analysis, aes(x = as.factor(RegionType), y = AverageSales, fill = as.factor(Group))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ ProductType + Style) +
  ggtitle("Sales by Region, Product, and Style for Each Group") +
  xlab("Region Type") +
  ylab("Average Sales") +
  scale_fill_discrete(name = "Group") +
  theme_minimal()

# 後驗概率與銷售額的相關性
data_with_group <- cbind(data, PosteriorMean = latent_probs$PosteriorMean, Group = latent_probs$Group)

correlation_analysis <- data_with_group %>%
  summarise(Correlation = cor(PosteriorMean, Sales))
print(correlation_analysis)

ggplot(data_with_group, aes(x = PosteriorMean, y = Sales)) +
  geom_point(aes(color = as.factor(Group)), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("Posterior Probability vs Sales") +
  xlab("Posterior Mean of z[i]") +
  ylab("Sales") +
  theme_minimal() +
  scale_color_discrete(name = "Group")


################
## 第五個模型 ##
################
jags_model <- "
model {
  for (i in 1:N) {
    Sales[i] ~ dt(mu[i], tau, df)  # t-distribution for robustness
    mu[i] <- beta_0 + beta_quarter[Quarter[i]] + beta_region[RegionType[i]] +
             beta_product[ProductType[i]] + beta_style[Style[i]] +
             beta_interaction[ProductType[i], Style[i], RegionType[i]]
  }

  # 熵的計算，用於模型擴展
  for (p in 1:2) {
    for (s in 1:3) {
      for (r in 1:2) {
        exp_beta_interaction[p, s, r] <- exp(beta_interaction[p, s, r])
      }
    }
  }

  for (p in 1:2) {
    for (s in 1:3) {
      for (r in 1:2) {
        prob[p, s, r] <- exp_beta_interaction[p, s, r] / sum(exp_beta_interaction[p, s, ])
        entropy[p, s, r] <- -prob[p, s, r] * log(prob[p, s, r] + 1e-10)
      }
    }
  }

  # Priors for fixed effects
  beta_0 ~ dnorm(0, 0.001)

  for (j in 1:4) {
    beta_quarter[j] ~ dnorm(mean_quarter[j], prec_quarter[j])
  }

  for (k in 1:2) {
    beta_region[k] ~ dnorm(mean_region[k], prec_region[k])
  }

  for (m in 1:2) {
    beta_product[m] ~ dnorm(mean_product[m], prec_product[m])
  }

  for (n in 1:3) {
    beta_style[n] ~ dnorm(mean_style[n], prec_style[n])
  }

  for (p in 1:2) {
    for (q in 1:3) {
      for (r in 1:2) {
        beta_interaction[p, q, r] ~ dnorm(0, 0.001)
      }
    }
  }

  # Variance priors
  tau <- 1 / (sigma^2)
  sigma ~ dunif(0, 100)
  df ~ dexp(1)
}
"

model <- jags.model(textConnection(jags_model), data = jags_data, inits = inits, n.chains = 3, n.adapt = 1000)
update(model, 10000)
samples <- coda.samples(model, variable.names = c("beta_0", "beta_quarter", "beta_region", "beta_product", "beta_style", "beta_interaction", "entropy"), n.iter = 5000)

summary(samples)

mcmc_samples <- do.call(rbind, lapply(samples, as.data.frame))
mcmc_samples <- as.data.frame(mcmc_samples)

# 款式效應
style_effects <- mcmc_samples %>%
  dplyr::select(starts_with("beta_style")) %>%
  pivot_longer(cols = everything(), names_to = "Style", values_to = "Effect")

# 品項定位效應
product_effects <- mcmc_samples %>%
  dplyr::select(starts_with("beta_product")) %>%
  pivot_longer(cols = everything(), names_to = "ProductType", values_to = "Effect")


ggplot(style_effects, aes(x = Effect, fill = Style)) +
  geom_density(alpha = 0.5) +
  ggtitle("Posterior Distributions of Style Effects by Region") +
  theme_minimal()

ggplot(product_effects, aes(x = Effect, fill = ProductType)) +
  geom_density(alpha = 0.5) +
  ggtitle("Posterior Distributions of Product Effects by Region") +
  theme_minimal()

gelman.diag(samples)
