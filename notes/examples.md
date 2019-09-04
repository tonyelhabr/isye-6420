---
title: 'Bayesian Theory and Examples'
author: ''
output:
  rmdformats::html_docco:
    keep_md: true
    highlight: kate
    toc: true
    toc_depth: 5
---





<style type="text/css">
body, td {
  font-size: 16px;
}
code.r {
  font-size: 12px;
}
pre {
  font-size: 12px;
}
ref {
  font-style: italic;
}
hide {
  display: none;
}
img.image-thumb {
    max-width: 100%;
}
</style>

## Theory

<i>
A note about notation: In general, we use capital letters to note a random variable (RV)
that has not yet taken on a value, and a lower-case letter to indicate an
"instatiation" or "manifestation" of the RV.
</i>


### Probability Distribution

These aren't all of the distributions, but these are brought up
in the examples that follow. (Note that $f(x)$ represents the [probability mass function](https://en.wikipedia.org/wiki/Probability_mass_function) (PMF)
for discrete distribution and the [probability density function](https://en.wikipedia.org/wiki/Probability_density_function) (PDF) for continuous
distribution

#### Binomial

$$
\begin{array}{c}
X \sim \mathcal{Bin}(n, p) \\
f(x)=\left(\begin{array}{l}{n} \\ 
{k}\end{array}\right) p^{X} q^{n-X}, k=0,1, \ldots, n \\
\operatorname{E}[X]=n p, \\
\operatorname{Var}(X)=n p q.
\end{array}
$$

This is used to model the number of successes $x$ in $n$ trials, 
where the probability of success is $p$. (It's a generalization of the Bernoulli distribution, which is limited to just 1 trial.)

#### Uniform

$$
\begin{array}{c}
X \sim \mathcal{U}(a, b) \\
f(x)=\frac{1}{b-a}, a \leq x \leq b, \\
\operatorname{E}[X]=\frac{a+b}{2}, \\
\operatorname{Var}(X)=\frac{(b-a)^{2}}{12}.
\end{array}
$$

This is commonly used if we only really know the upper and lower bounds of a distribution, and not much else.

#### Normal

$$
\begin{array}{c}
X \sim \mathcal{N}(\mu, \sigma^2) \\
f(x)=\frac{1}{\sqrt{2 \pi \sigma^{2}}} \exp \left[\frac{-(x-\mu)^{2}}{2 \sigma^{2}}\right] \\
\operatorname{E}[X]=\mu, \\
\operatorname{Var}(X)=\sigma^{2}.
\end{array}
$$

This is the "mother" of all distributions.
The sum of the values coming from any distribution will eventually 
(approximately) follow the normal distribution.
(See the Central Limit Theorem.)

<hide>
#### Gamma

$$
\begin{array}{c}
X \sim \mathcal{Ga}(\alpha, \beta) \\
f(x)=\frac{\beta^{\alpha} x^{\alpha-1} e^{-\beta x}}{\Gamma(\alpha)}, x \geq 0 \\
\operatorname{E}[X]=\frac{\alpha}{\beta}, \\
\operatorname{Var}(X)=\frac{\alpha}{\beta^{2}}.
\end{array}
$$

Note that $\Gamma(\alpha) = \int_{0}^{\infty} t^{\alpha-1} e^{-t} dt.$

The Gamma distribution is a generalization of two other well-known distributions---the exponential and chi-squared distributions (which themeselves
draw upon attributes of other distributions). It shows up in a couple of conjugate pairs (to be discussed).
</i>

#### Beta

$$
\begin{array}{c}
X \sim \mathcal{Be}(\alpha, \beta) \\
f(x)=\frac{1}{\operatorname{B}(\alpha,\beta)} x^{\alpha - 1} (1-x)^{\beta-1}, 0 \leq x \leq 1, \alpha, \beta > 0. \\
\operatorname{E}[X]= \frac{\alpha}{\alpha + \beta}, \\
\operatorname{Var}(X)=\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}
\end{array}
$$

Note that $\operatorname{B}(\alpha,\beta) = \frac{\Gamma(\alpha) \Gamma(\beta)}{\Gamma(\alpha + \beta)}$ and that $\Gamma(\alpha) = \int_{0}^{\infty} t^{\alpha-1} e^{-t} dt.$

The Beta distribution is a generalization of the Gamma distribution, which
is itself a generalization of two other well-known distributions---the exponential and chi-squared distributions (which themeselves
draw upon attributes of other distributions).
The Beta probably the best distribution for modeling other distributions.

### Bayes' Theorem

<ref>
"Bayesian probability is an interpretation of the concept of probability, in which ... probability is interpreted ... as quantification of a personal belief." [^1]
</ref>

$$
\begin{array}{ccc}
\text{ Posterior Probability } & = & \frac{(\text{ Likelihood })(\text{ Prior Probability })}{\text{ Marginal Likelihood } } \\
P(H|E) & = & \frac{P(E|H) P(H)}{P(E)}
\end{array}
$$

<i>
Note the notation: We use $H$ to indicate "hypothesis" and $E$ to indicate evidence. Sometimes evidence is also referred to as "data" ($D$) in the literature.
</i>

Also, note that the Bayes' formula is really just a different way of stating the
conditional probability of events. (i.e. The probability of $A$ given $B$ is $P(A|B) = \frac{P(B|A) P(B)}{P(A)}$.)

Regarding the components of Bayes' formula: 

+ The prior $P(H)$ is the probability of $H$ before observing the data. Referring back to the definition of Bayes probability, it is the "personal" component. It is "elicited" by the person observing the data/performing the experiment, and its reliability is dependent on the knowledge of this person.

+ The likelihood $P(E|H)$ is the evidence about $H$ provided by the data.

+ The marginal $P(E)$ (i.e. the "normalizing constant") is the probability of the observing the data, accounting for all possibile hypotheses. It is usually what causes difficulty in calculations. We can avoid
difficulty in calculations if our problem can be modeled with conjugate pairs
(to be discussed). Also, note that the marginal may also be referred to as the
 predictive prior. A more in-depth discussion of this discussion would clarify
 why this is.

+ The posterior $P(H|E)$ is the probability that $H$ is true after the data is considered.

### Bayesian vs. Frequentist/Classical [^2]

<ref>
"Disagreements are in the nature of (1) model parameters and (2) use of conditioning:
</ref>

+ <ref>Frequentists: (1) Parameters are <b>fixed numbers</b>. (2) Inference involves <b>optimization.</b></ref>

+ <ref>Bayesians: Parameters are <b>random variables</b>. (2) Inference involves <b>integration</b>."</ref>

(See the Beta-Binomial example below.)

### Conjugacy

<ref>
"Conjugacy occurs when the posterior distribution is in the same family of probability density functions as the prior belief, but with new parameter values, which have been updated to reflect what we have learned from the data." [^3]
</ref>

#### Beta-Binomial

<ref>
"Suppose we perform an experiment and estimate that the data comes from a binomial distribution $\mathcal{Bin}(n, p)$ with <b>known</b> $n$ (number of trials) and <b>unknown</b> $p$ (probability). (This is the likelihood.)
Let's say that we have a prior belief that $p$ can be modeled with the Beta distribution $\mathcal{Be}(\alpha, \beta)$ (with $\alpha, \beta$ that we choose, presumably with our "expertise"). (This is the prior.)
In the experiment we observed $x$ successes in $n$ trials. Then Bayes' rule implies that our new belief about the probability density of $p$---the posterior distribution, of $p | x$---is also the Beta distribution, with "updated" parameters." [^3]
</ref>

$$
\begin{array}{c}
p | x \sim \mathcal{Be}(\alpha + x, \beta + n - x).
\end{array}
$$

We can formalize this problem set-up as follows.

$$
\begin{array}{rclcl}
\text{ Likelihood } & : & x | p & \sim & \mathcal{Bin}(n, p) \\
\text{ Prior } & : & p & \sim & \mathcal{Be}(\alpha, \beta) \\
\text{ Posterior } & : & p | x & \sim & \mathcal{Be}(\alpha + x, \beta + n - x)
\end{array}
$$

<i>
Note that it is common to "simplify" notation, e.g. $P(p|x)$ is express simply as $p|x$.
</i>

The posterior mean reflects an "update" to the prior mean 
given $x$ and $n$.

$$
\begin{array}{rcl}
\operatorname{E}[X] & \rightarrow & \operatorname{E}[p|X] \\
\frac{\alpha}{\alpha + \beta} & \rightarrow & \frac{\alpha + x}{\alpha + \beta + n}.
\end{array}
$$

##### Example

We want to estimate the probability $p$ that a coin falls heads up. After $n = 10$ flips, we observe $X = 0$ heads.
we can model the likelihood as a $\mathcal{Bin}(n = 10, p)$ (where $p$ is unknown).
A "flat" prior for this kind of experiment is a $\mathcal{U}(0, 1)$ distribution, which happens to be equivalent to $\mathcal{Be}(1, 1)$.
What does the Beta-Binomial conjugate pair lead us to conclude about $p$?

Using $\alpha=1, \beta=1$ for our prior distribution, the posterior probability is  $\mathcal{Be}(1 + (0), 1 + (10) - (0)) = \mathcal{Be}(1, 11)$. Then the posterior mean (i.e. average or expectation) is $\hat{p} = \frac{(1)}{(1)+(11)} = \frac{1}{12}$ (because the expectation of the $\mathcal{Be}$ distribution is $\operatorname{E}[X] = \frac{\alpha}{\alpha+\beta}$.) Note that 

If we had taken the Frequentist approach (which does not incorporate priors), then
we would have concluded that $\hat{p} = \frac{X}{n} = \frac{0}{10} = 0$.

Note that the Bayes estimation of the posterior mean would be much closer to $p = 0.5$ if we had used a stronger prior (which is more realistic for something like
coin flips, where it is physically very difficult to create a "rigged" coin). For example, if we had used the prior $\mathcal{Be}(1000, 1000)$, then our posterior estimate of the mean
would have been $\hat{p} = \frac{(10000)}{(10000) + (11)} \approx 0.5$.

#### Normal-Normal

For the Normal-Normal conjugate pair, we assume that the data comes 
from a normal distribution with <b>known</b> variance $\sigma^2$ and <b>unknown</b> mean $\mu$, which we want to estimate. Also, we must "elicit"
value for mean $\mu_0$ and variance $\sigma_0^2$ for the prior distribution (which is
itself normal). (This is like how we choose $\alpha$ and $\beta$ for the prior for the Beta-Binomial conjugate pair.)

$$
\begin{array}{rclcl}
\text{ Likelihood } & : & x | \mu & \sim & \mathcal{N}(\mu, \sigma^2) \\
\text{ Prior } & : & \mu & \sim & \mathcal{N}(\mu_0, \sigma_0^2) \\
\text{ Posterior } & : & \mu | x & \sim & \mathcal{N}(\frac{\mu_0 \sigma^2 + n x \sigma_0^2}{\sigma^2 + n \sigma_0^2}, \frac{\sigma^2 \sigma_0^2}{\sigma^2 + n \sigma_0^2})
\end{array}
$$

If you look closely, you'll see that the posterior mean $\mu$ is actually
a weighted average of the prior mean $\mu_0$ and the observed mean $\frac{x}{n}$.

##### Example

Joe models his IQ as $X \sim \mathcal{N}(\mu, 80)$.
The distribution of IQs of students at Joe's university is
$\mathcal{N}(110, 120)$.
Joe takes a single IQ test and scores $98$.


$$
\begin{array}{rclcl}
\text{ Likelihood } & : & x | \mu & \sim & \mathcal{N}(\mu, 80) \\
\text{ Prior } & : & \mu & \sim & \mathcal{N}(110, 120) \\
\text{ Posterior } & : & \mu | x & \sim & \mathcal{N}(\frac{(110) (80) + (1) (98) (120)}{(80) + (1) (120)}, \frac{(80) (120)}{(80 + (1) (120)}) \approx \mathcal{N}(102.8, 48)
\end{array}
$$

### General Bayes Computation

So we have seen how posterior distributions for parameters can be generated via conjugacy.
Unfortunately, very few settings can be modeled "well" with 
conjugate pairs--which are nice because they have analytical, closed form solutions. 
Probably the first "computational" approach to try is a grid approximation. 
(In fact, this approach is useful for verifying closed-form problems.)
However, this approach tends to be too simplistic and is usually not very "useful". 
The most popular computation approach--and one that is used a lot in practice--is MCMC.
This is where "specialized" Bayesian software (e.g. OpenBUGS, JAGS, Stan) come into
play. However, these are not "native" to `R`. (There are `R` packages that
"wrap" the functionality of these software, but I wouldn't consider this "native".)

If you want to be as sophisticated
as possible without stepping outside of the realm or `R`, then the technique
to go with is is Laplace approximation, which is also known as quadratic approximation.
By providing a means to approximate parameters for a "typical" linear regression,
Laplace approximation can be used as the "engine" for "Bayesian linear regression".


#### Laplace Approximation [^4]

<ref>
"I’ll demonstrate the use of quadratic approximation with a model for an individual team’s score differential. The first step is to gather the score differential for each team in each season...
for NFL seasons from 2009 to 2017.
</ref>



```r
# Load the tidyverse
# install.packages("tidyverse")
library(tidyverse)

# Use the nflWAR package to get summaries of team performances for each season,
# first install devtools if you don't have then install nflWAR:
# install.packages("devtools")
# devtools::install_github("ryurko/nflWAR")
library(nflWAR)

# Use the get_season_summary() function to create a dataset that has a row
# for team-season with a column for their score differential from 2009 to 2017:
team_summary_df <- get_season_summary(2009:2017)
```

<ref>
Before going into the regression example with a predictor, it’s worthwhile to first demonstrate quadratic approximation by just modeling the score differential with a Gaussian.
</ref>


```r
# Create a histogram displaying the score differential distribution using
# the the density instead so it's consistent for later: (also why I'm assigning
# it to a variable name, so we can add our approximated posterior layer later)
score_diff_hist <- team_summary_df %>%
  ggplot(aes(x = Total_Score_Diff)) +
  geom_histogram(aes(y = ..density..), bins = 20) + theme_bw() + 
  scale_x_continuous(limits = c(-300, 300)) +
  # Always label:
  labs(x = "Score differential", y = "Density",
       title = "Distribution of individual team score differential in each season from 2009-17 ",
       caption = "Data accessed with nflscrapR") +
  # Just some cosmetic settings: (really should do a custom theme at the 
  # beginning instead of repeating this)
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14))

# Display the plot
score_diff_hist
```

<img src="examples_files/figure-html/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

<ref>
"From this histogram we have support for using a Gaussian, since it’s roughly symmetric and unimodal.<br/><br/>
We assume then that the general model for an observed score differential $S_i$, where $i = 1, \dots, 288$
is the index for an observed team-season score differential, follows a Gaussian distribution:
</ref>

$$
\begin{array}{c}
S_i \sim \mathcal{N}(\mu, \sigma^2)
\end{array}
$$

<ref>
"Because we’re approaching this problem from a Bayesian perspective, we really want to consider an infinite number of possible Guassians based on all of the possible combinations for our parameters $\mu$ and $\sigma$. Of course we’re not going to actually look at every possible value for these parameters, Laplace approximation will accomplish what we need. But, ... we’re interested in the entire posterior distribution - which in this case is a posterior distribution of distributions. So for a full model of score differential, we need priors for both  $\mu$ and $\sigma$.
</ref>

$$
\begin{array}{c}
S_i \sim \mathcal{N}(\mu, \sigma^2) \\
\mu \sim \mathcal{N}(0, 100^2) \\
\sigma \sim \mathcal{U}(0,200)
\end{array}
$$
<ref>
"This is a really naive model, just to demonstrate Laplace approximation. Next we’ll implement the [some] code ... which uses the `optim()` function in R to find the mode of the posterior and also returns the Hessian matrix for the curvature (following the observed information definition above). For ease, we’ll ... [use] the `{mvtnorm}` package for sampling from the approximated posterior. (This sampling isn’t really necessary but it’s good practice to start sampling from the posterior). We’re working with two parameters here, so the resulting posterior approximation is a multivariate Gaussian, hence why we use the `rmvnorm()` function since we’re not only working with each parameter’s mean and variance, but their covariance as well. (In this example the covariance is essentially 0, but that won’t typically be the case).
</ref>


```r
# Function to perform simple Laplace approximation
# @param log_posterior_n Function that returns the log posterior given a vector
# of parameter values and observed data.
# @param init_points Vector of starting values for parameters.
# @param n_samples Number of samples to draw from resulting approximated 
# posterior distribution for then visualizing.
# @param ... Any additional arguments passed to the log_posterior_fun during 
# the optimization.
# @return Samples from approximated posterior distribution.
laplace_approx <- function(log_posterior_fn, init_points, n_samples, ...) {
  
  # Use the optim function with the input starting points and function to
  # evaluate the log posterior to find the mode and curvature with hessian = TRUE
  # (note that using fnscale = -1 means to maximize). We use the ... here to 
  # allow for flexible name of the data in log_posterior_fn:
  fit <- optim(init_points, log_posterior_fn, 
               # Use the same quasi-Newton optimization method as McElreath's
               # map function in his rethinking package:
               method = "BFGS",
               # using fnscale = -1 means to maximize
               control = list(fnscale = -1), 
               # need the hessian for the curvature
               hessian = TRUE, 
               # additional arguments for the function we're optimizing
               ...)
  
  # Store the mean values for the parameters and curvature to then generate
  # samples from the approximated normal posterior given the number of samples
  # for both of the parameters, returning as a data frame:
  param_mean <- fit$par
  # inverse of the negative hessian for the covariance - same as the use 
  # of the observed information in the intro
  param_cov_mat <- solve(-fit$hessian)
  # sample from the resulting joint posterior to get posterior distributions
  # for each parameter:
  mvtnorm::rmvnorm(n_samples, param_mean, param_cov_mat) %>%
    data.frame()
}

# Function to calculate log posterior for simple score differential model
# @param params Vector of parameter values that are named "mu" and "sigma".
# @param score_diff_values Vector of observed score differential values
score_diff_model <- function(params, score_diff_values) {
  
  # Log likelihood:
  sum(dnorm(score_diff_values, params["mu"], params["sigma"], log = TRUE)) +
    # plus the log priors results in log posterior:
    dnorm(params["mu"], 0, 100, log = TRUE) + 
    dunif(params["sigma"], 0, 200, log = TRUE)

}
```

<ref>
"With these functions we’ll now approximate the posterior, considering very naive starting values for our optimization with $\mu = 200$ and $\sigma = 10$. (This is just to demonstrate, I have no reason for choosing these values).
</ref>


```r
init_samples <- laplace_approx(log_posterior_fn = score_diff_model,
                               init_points = c(mu = 200, sigma = 10),
                               n_samples = 10000,
                               # Specify the data for the optimization!
                               score_diff_values = team_summary_df$Total_Score_Diff)
```


<ref>
"We can easily view the joint distribution of our posterior samples with their respective marginals on the side.
</ref>


```r
# install.packages("cowplot")
library(cowplot)

# First the joint distribution plot with points for the values:
joint_plot <- ggplot(init_samples, aes(x = mu, y = sigma)) + 
  geom_point(alpha = .3, color = "darkblue")  +
  labs(title = "Posterior distribution for model parameters with marginals") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16))

# Marginal density for mu along x-axis using the cowplot function axis_canvas:
mu_dens <- axis_canvas(joint_plot, axis = "x") +
  geom_density(data = init_samples, aes(x = mu), fill = "darkblue",
               alpha = 0.5, size = .2)

# Same thing for sigma but along y and with coord_flip = TRUE to make it vertical:
sigma_dens <- axis_canvas(joint_plot, axis = "y", coord_flip = TRUE) +
  geom_density(data = init_samples, aes(x = sigma), fill = "darkblue",
               alpha = 0.5, size = .2) +
  coord_flip()

# Now generate by adding these objects to the main plot:
# Need grid:
# install.packages("grid")
joint_plot_1 <- insert_xaxis_grob(joint_plot, mu_dens, 
                                  grid::unit(.2, "null"), position = "top")
joint_plot_2 <- insert_yaxis_grob(joint_plot_1, sigma_dens, 
                                  grid::unit(.2, "null"), position = "right")
ggdraw(joint_plot_2)
```

<img src="examples_files/figure-html/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />


<ref>
"From this we can clearly see the Gaussian approximations for both of our parameters with the distribution of 
$\mu$ centered around 0 and the distribution for $\sigma$ centered around 100. Remember this represents the distributions for the parameters of our score differential model. So sampling from this posterior means we are generating different Gaussian distributions for the score differential. ... [T]he code below uses these 10000 values from `init_samples()` for each parameter, and then samples 10000 values from distributions using these combinations of values to give us our approximate score differential distribution.
</ref>
 

```r
# R is vectorized so can just give it the vector of values for the Gaussian
# distribution parameters to create the simulated score-diff distribution:
sim_score_diff_data <- data.frame(Total_Score_Diff = 
                                    rnorm(10000, mean = init_samples$mu,
                                          sd = init_samples$sigma))

# Create the histogram from before but now add this posterior density on top:
score_diff_hist + geom_density(data = sim_score_diff_data,
                               aes(x = Total_Score_Diff), 
                               fill = NA, color = "darkblue")
```

<img src="examples_files/figure-html/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

<ref>
"Not bad, but then again this was a pretty easy example. You could also use the grid approximation from the last post to accomplish the same task but it starts to be annoying to work with even with just two parameters, quickly becoming unbearable with more and more parameters.
</ref>


#### Bayesian Linear Regression with Laplace Approximation [^4]

<ref>
"Ok, so we’ve seen how easy it is to use Laplace approximation for modeling score differential… but that didn’t tell us anything interesting. What we’re really interested in, is the relationship between some predictor variable(s) and a team’s score differential. For instance, how is a team’s passing performance related to their score differential? What about their pass defense? Is this more or less related to score differential than rushing?
</ref>

<ref>
"... [W]e’re now going use both efficiency differentials as predictors. Remember, all we need to fully describe a Gaussian distribution is the mean and variance. So in order to incorporate these into the our model we’ll define the mean, $\mu$ as a function of both predictors explicitly. For simplicity we’ll return to using $i$ to denote a single combination of team $t$ and season $y$ for the 286 team-season combinations that we have, and we’ll also let $P_i$ and $R_i$ be representations of passing and running efficiency:
</ref>

$$
\begin{array}{c}
S_i \sim \text{N}(\mu_i, \sigma^2) \\
\mu_i = \alpha + \beta_{P}P_i + \beta_R R_i \\
\alpha \sim \text{N}(0, 100^2) \\
\beta_P \sim \text{N}(0, 100^2) \\
\beta_R \sim \text{N}(0, 100^2) \\
\sigma \sim \text{Uniform}(0,200)
\end{array}
$$


<ref>
"Now the mean of the score differential distribution _depends on the predictors_
as denoted with the subscript $i$. By denoting the relationship between $\mu_i$ 
and each of the predictors with an $=$ sign, we're saying that $\mu_i$ is no 
longer a parameter, its just a function of our actual parameters of interest.
The parameters we really care about here are $\alpha$, $\beta_P$, and $\beta_R$.
Of course all of this should look familiar if you've seen [linear regression](https://en.wikipedia.org/wiki/Linear_regression)
before. The $\alpha$ is our model's intercept, representing the _expected_ score
differential for when both the passing and rushing efficiency differentials 
are 0. And each of the $\beta$s represent the change in _expected_ score 
differential when the pass/rush efficiency differential increases by 1 unit, after
accounting for the other type of efficiency differential (rush/pass). So each of
these are given priors (along with $\sigma$ again), and we don't need to specify 
a prior for $\mu_i$ since it is explained completely by these parameters. These
priors for now are pretty weak, but both $\beta$ priors centered at 0 are just
saying, "I don't know if there's a relationship between these effiency 
differentials and score differential, good or bad are equally likely.
</ref>

<ref>
"To actually approximate this model..."
</ref>

--------------------------------------------------------------------------------

[^1]: https://en.wikipedia.org/wiki/Bayesian_probability.

[^2]: Tony's Bayes Statistics class.

[^3]: https://statswithr.github.io/book/bayesian-inference.html#conjugacy

[^4]: http://www.stat.cmu.edu/~ryurko/post/bayesian-baby-steps-normal-next-steps

