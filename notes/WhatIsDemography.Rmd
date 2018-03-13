# 

## Intro

### How does UW breakdown there Demography group  
 - Pop Health  
 - Meathods and Measurement  
 - Migration and Mobility  
 - Pop and Environment  
 - Child/Family Well Being  

## Balancing Equations  
 - $K_{2016} = K_{2015} + B_{2105} - D_{2015} + \epsilon + \zeta$  
 - Geometric model  
   - $K_{2016} = K_{2015} (1 + \frac{B_{2015}}{K_{2015}} - \frac{D_{2015}}{K_{2015}})  
   - $K_{2016} = K_{2015} (1 + R_{2015})$  
   - $\frac{K_t}{K_0} = (1 + R)^t$  
   - $R = (\frac{k_t}{K_0})^{1/t} - 1$  
   - $K_t = K_0 (1 + \frac{R}{m})^{mt}$  
   - $lim_{m \rightarrow \infty} (1 + \frac{R}{m})^m = e^R$  
   - $\therefore  K_t = K_0 e^{Rt}$  
   - $\text{log}(K_t) - \text{log}(K_0) = Rt$  
   - More Generally  
   - $R = \frac{\text{log}(K_{t_2}) - \text{log}(K_{t_1})}{t_2 - t_1}$  
     - Looks like a derivative!  
 - Another perspective  
   - $\frac{dK}{dt} = RK$  
   - $t^{\text{double}} = \frac{log(2)}{R}$  

```{R}
K <- c("1950"=2558000000, "2000"=6089000000) # population at two time points  
R <-  diff(log(K)) / (2000-1950) # calculate R from given equations above  
log(2) / R # calculate the doubling time with R
```

 - Logistic model
   - $\frac{dK}{dt} = RK(1 - \frac{K}{C})$  
   - Where $C$ is the carrying capacity  
   - $K_t = \frac{C}{1 + Ae^{-Rt}}$  
   - Where $ A = \frac{C-K_0}{K_0}$  

 - Question???  
   - Given that you have 100 fish and a carrying capcity of 500 fish how long till 250 fish?

```{R}
K0 <- 100
Kt <- 250
C <- 500
A <- (C - K0) / K0
# 250 = 500 / (1 + A * exp(-.1t))
# 2 = 1 + A * exp(-.1t)
# A^-1 = exp(-.1t)
# log(A^-1) = -.1t
# log(A^-1) / -.1 = t
log(A^-1) / -.1
```

## Period Perspective  
 - Period Person Years Lived = PPYL = PY  
   - The amount of time (years) that individuals in the population spend alive between two times  
 - Stable (age composition is not changing) vs Stationary (R=0)  
 - Cohort Person Years Lived = $B e_0 T = K b e_0 T$  
 - Period Person Years lived = $K T$  
 - If $CPYL == PPYL$ then $b e_0 = 1$ ie the stationary pop identity  
 - If $b$ and $e_0$ are both changing then you can have a stationary population that is not stable  

