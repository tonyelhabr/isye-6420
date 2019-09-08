
library(R2OpenBUGS)
model <- function() {
  E ~ dbern(0.5)
  c <- max(E, 0.5)

  works ~ dbern(c)
}

t = 0.5

data <- list(p = c(exp(-t), exp(-2*t), exp(-t / 2), exp(-t/3), exp(-t)),
             works = 1)

inits <- NULL

out <- bugs(data = data,
            inits = inits,
            parameters.to.save = c("E", "c"),
            model.file = model,
            digits = 5,
            n.chains = 3,
            n.burnin = 100,
            n.iter = 1000,
            # OpenBUGS.pgm=OpenBUGS.pgm,
            # WINE = WINE,
            # WINEPATH = WINEPATH,
            # useWINE=T,
            # debug = T
            debug = F
)

out
print(out$mean$E)
print(out$sd$E)
