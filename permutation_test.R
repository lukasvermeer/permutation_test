# Permutation Test Test
#
# A test to test the permutation test by testing many permutations of the permutation test.
#
# This script compiles the C code, generates some random data, runs the test, and plots histograms.
#
# You need GSL. Installed easily on a Mac with Brew (http://brew.sh/).
#

# WHERE ARE TEH FILEZ?!
p <- "/Users/lvermeer/Documents/GitHub/permutation_test/"

# PLZ WAIT COMPILING!!1
system2("gcc", args=paste("-DDEBUG=0 -O2 -Wall ",p,"permutation_test.c -o ",p,"permutation_test -lgsl -lgslcblas -lpthread -lm",sep=""))
system2("chmod", args=paste("+x ",p,sep=""))
  
p_dist <- c()
# I WILL  ^^^ STORE MAH DATA HERE

for (i in 1:100) {
  # IM JUST MAKING THIS SHIT UP AS I GO ALONG
  a <- rnorm(100, mean = 0)
  b <- rnorm(100, mean = 0)

  # I CAN HAZ P VALUES?
  p_vals <- as.numeric(system2(paste(p,"permutation_test",sep=""), input = c(paste(a, collapse = " "), paste(b, collapse = " ")), stdout = T))
  p_dist <- append(p_dist, mean(p_vals))
}

# FANZY PLOTZ
par(mfrow=c(2,1))

hist(p_vals, breaks=20, main = "Distribution of the last random sample", xlab = "p-value\nt-test (red) vs. average permutation test value (black)")
t <- t.test(a,b)
abline(v=t$p.value, col="red", lwd=2)
abline(v=mean(p_vals), col="black", lwd=2)

hist(p_dist, xlim=c(0,1), main = "Distribution of p-values for many random samples", xlab = "p-value")