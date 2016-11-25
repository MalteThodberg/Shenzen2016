library(knitr)
library(tidyverse)
library(forcats)
library(broom)

set.seed(2016)
options(digits=3)

#' Linear regression
d <- tibble(x=1:20,
						y=x+rnorm(20))

ggplot(d, aes(x=x, y=y)) +
	geom_point(alpha=0.75) +
	geom_smooth(method="lm", se=FALSE, alpha=0.75, color="red") +
	theme_minimal()

#' Setup Example simulations

# Simple group experiments
twogroup <- tibble(control=rnorm(10, 5, 1),
											treated=rnorm(10, 7, 1)) %>%
	gather(key="Group", value="Expression")

threegroup <- tibble(group1=rnorm(10, 4, 1),
									 group2=rnorm(10, 9, 1),
									 group3=rnorm(10, 6, 1)) %>%
	gather(key="Group", value="Expression")

# Interactions
noint <- tibble(con_wt=rnorm(10, 4, 1),
								con_mut=rnorm(10, 7, 1),
								trt_wt=rnorm(10, 6, 1),
								trt_mut=rnorm(10, 9, 1)) %>%
	gather(key="Group", value="Expression") %>%
	separate(Group, into=c("Treatment", "Condition"), remove=FALSE) %>%
	mutate(Condition=fct_relevel(Condition, "wt"),
				 Group=fct_relevel(Group, c("con_wt", "con_mut", "trt_wt", "trt_mut")))

withint <- tibble(con_wt=rnorm(10, 4, 1),
								con_mut=rnorm(10, 7, 1),
								trt_wt=rnorm(10, 6, 1),
								trt_mut=rnorm(10, 12, 1)) %>%
	gather(key="Group", value="Expression") %>%
	separate(Group, into=c("Treatment", "Condition"), remove=FALSE) %>%
	mutate(Condition=fct_relevel(Condition, "wt"),
				 Group=fct_relevel(Group, c("con_wt", "con_mut", "trt_wt", "trt_mut")))

# Batches
batch <- tibble(batchA_group1=rnorm(10, 4, 1),
								batchA_group2=rnorm(10, 6, 1),
								batchB_group1=rnorm(10, 12, 1),
								batchB_group2=rnorm(10, 14, 1),
								batchC_group1=rnorm(10, 8, 1),
								batchC_group2=rnorm(10, 10, 1)) %>%
	gather(key="tmp", value="Expression") %>%
separate(tmp, into=c("Batch", "Group"))


bad <- tibble(group1=rnorm(10, 4, 1),
							group2=rnorm(10, 12, 1)) %>%
	gather(key="Group", value="Expression") %>%
	mutate(Batch=ifelse(Group=="group1", "BatchA", "BatchB"))

#' Plots and linear model fits

# First plot as template
p <- ggplot(twogroup, aes(x=Group, y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed") +
	scale_color_brewer(palette="Set1") +
	theme_minimal()

# Simple groups
kable(tidy(summary(lm(Expression~Group, data=twogroup))))
p
kable(tidy(summary(lm(Expression~Group, data=threegroup))))
p %+% threegroup

# Interactions
kable(tidy(summary(lm(Expression~Treatment*Condition, data=noint))))
p %+% noint + scale_color_manual(values=c("grey", "red", "blue", "purple"))
kable(tidy(summary(lm(Expression~Treatment*Condition, data=withint))))
p %+% withint + scale_color_manual(values=c("grey", "red", "blue", "purple"))

# batch
kable(tidy(summary(lm(Expression~Batch+Group, data=batch))))
p %+% batch + facet_wrap(~Batch)
kable(tidy(summary(lm(Expression~Batch+Group, data=bad))))
p %+% bad + facet_wrap(~Batch)

#' Use QR decomposition to check for full rank

# Extract model matrix
batch_mod <- model.matrix(Expression~Batch+Group, data=batch)
bad_mod <- model.matrix(Expression~Batch+Group, data=bad)

# Number of columns
ncol(batch_mod)
ncol(bad_mod)

# Independent columns
qr(batch_mod)$rank
qr(bad_mod)$rank

# Indendent columns must be larger than total number of columns
qr(batch_mod)$rank >= ncol(batch_mod)
qr(bad_mod)$rank >= ncol(bad_mod)

#' Simulate data directly from the model matrix

simulate_lm <- function(mod, Beta){
	stopifnot(ncol(mod) == length(Beta))
	mod %*% Beta + rnorm(ncol(mod))
}

simulate_lm(mod=model.matrix(~Treatment*Condition, data=withint),
						Beta=c(1,2,3,4)) %>%
	head
