# FULL SYNTHETIC fixed-protocol evaluation; base R only. No package installation.
# The bounded foreground Python launcher records commands, seeds and hashes.
options(warn=1)
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
full_design <- list(n=c(2048L,8192L,32768L,131072L),gamma=c(.25,.5,1),
 R=2000L,B=19L,tau=.5,trim=.1,alpha=.05,seed=2026111907L)
smoke_design <- list(n=c(2048L,131072L),gamma=.5,R=2L,B=19L,
 tau=.5,trim=.1,alpha=.05,seed=2026111908L)
get_design <- function(profile) {
 if(profile=="full") full_design else if(profile=="smoke") smoke_design else stop("Unknown profile")
}

clayton_one <- function(n) {
  # Gamma(shape=1,rate=1) frailty is Exp(1); exact Clayton(theta=1).
  w <- rexp(n)
  u <- w / (w + rexp(n))
  v <- w / (w + rexp(n))
  cbind(u = u, v = v)
}

make_sample <- function(n, gamma, tau) {
  location_from_uniforms(clayton_one(n),gamma,tau)
}

location_from_uniforms <- function(uv,gamma,tau) {
  n <- nrow(uv)
  k <- as.integer(floor(sqrt(n)))
  delta <- gamma / sqrt(k)
  stopifnot(delta > 0, delta < .5)
  b <- as.integer(floor(tau * n))
  x <- uv[, "u"]
  ii <- seq.int(b + 1L, n)
  x[ii] <- x[ii] - delta
  list(x = x, y = uv[, "v"], k = k, Delta = delta, b = b)
}

trim_indices <- function(n, trim) {
  seq.int(as.integer(ceiling(trim * n)), as.integer(floor((1-trim)*n)))
}

cusum_date <- function(x, trim) {
  n <- length(x)
  cs <- cumsum(x)
  jj <- trim_indices(n, trim)
  objective <- abs(cs[jj] - jj * cs[n] / n)
  jj[which.max(objective)]  # which.max picks the earliest exact tie.
}

support_date <- function(x, trim) {
  n <- length(x)
  neg <- which(x < 0)
  fallback <- length(neg) == 0L
  raw <- if (fallback) floor(n/2) else neg[1L] - 1L
  limits <- range(trim_indices(n, trim))
  j <- as.integer(max(limits[1L], min(limits[2L], raw)))
  list(j = j, clipped = j != raw, fallback = fallback)
}

joint_indicators <- function(x, y, j, k) {
  n <- length(x)
  stopifnot(length(y) == n, j >= 1L, j < n, !anyNA(x), !anyNA(y))
  ii <- seq_len(j)
  jj <- seq.int(j + 1L, n)
  rx <- numeric(n)
  rx[ii] <- rank(x[ii], ties.method = "first") / (j + 1)
  rx[jj] <- rank(x[jj], ties.method = "first") / (n - j + 1)
  ry <- rank(y, ties.method = "first") / (n + 1)
  as.integer(rx <= k/n & ry <= k/n)
}

integral_statistic <- function(ind) {
  # On [i/n,(i+1)/n), J(s)=J_i. Integrate (J_i-s*A)^2 exactly.
  n <- length(ind)
  a <- sum(ind)
  if (a == 0) return(0)
  ji <- c(0, cumsum(ind))[seq_len(n)]
  mid <- (seq_len(n) - 0.5)/n
  (mean((ji - mid*a)^2) + a^2/(12*n^2)) / a
}

drift_metrics <- function(i_est, i_oracle, j, b, k, Delta) {
  n <- length(i_est)
  tau <- b/n
  s <- (0:n)/n
  q <- pmin(s/tau, 1) - s
  dcount <- c(0, cumsum(i_est - i_oracle))
  total <- dcount[n+1L]
  amplitude <- 0.25 * max(j-b, 0) * Delta / sqrt(k)
  right <- (dcount - s*total)/sqrt(k) + amplitude*q
  left <- (dcount[seq_len(n)] - s[-1L]*total)/sqrt(k) + amplitude*q[-1L]
  mid <- (seq_len(n)-0.5)/n
  qmid <- pmin(mid/tau, 1)-mid
  rmid <- (dcount[seq_len(n)] - mid*total)/sqrt(k) + amplitude*qmid
  slope <- -total/sqrt(k) + amplitude*ifelse(mid < tau, 1/tau-1, -1)
  l2sq <- mean(rmid^2) + sum(slope^2)/(12*n^3)
  c(L2_residual = sqrt(l2sq), sup_residual = max(abs(c(right, left))),
    predicted_amplitude = amplitude)
}

write_append <- function(x, path) {
  for (nm in names(x)) if (is.double(x[[nm]])) x[[nm]] <- sprintf("%.17g",x[[nm]])
  write.table(x, path, sep = ",", row.names = FALSE, col.names = !file.exists(path),
              append = file.exists(path), quote = TRUE)
}

seed_row <- function(cell, replication, role, state) {
  out <- data.frame(cell = cell, replication = replication, role = role)
  for (i in seq_along(state)) out[[paste0("seed", i)]] <- state[i]
  out
}

cp_interval <- function(x, n, error) {
  if (x == 0L) lo <- 0 else lo <- qbeta(error/2, x, n-x+1)
  if (x == n) hi <- 1 else hi <- qbeta(1-error/2, x+1, n-x)
  c(lo, hi)
}

summarize_results <- function(dat, cfg, outdir) {
  cells <- unique(dat$cell)
  methods <- c("oracle", "cusum", "support")
  rej_rows <- list(); diff_rows <- list(); secondary_rows <- list(); frequency_rows <- list()
  for (cell in cells) {
    dd <- dat[dat$cell == cell, ]
    nr <- nrow(dd)
    for (method in methods) {
      yy <- dd[[paste0("reject_", method)]]
      ph <- mean(yy)
      point <- cp_interval(sum(yy), nr, .05)
      simultaneous <- cp_interval(sum(yy), nr, .025/(length(cells)*3))
      rej_rows[[length(rej_rows)+1L]] <- data.frame(
        cell = cell, n = dd$n[1L], gamma = dd$gamma[1L], method = method,
        trials = nr, rejections = sum(yy), estimate = ph,
        MCSE = sqrt(ph*(1-ph)/nr), CI95_low = point[1L], CI95_high = point[2L],
        simultaneous_low = simultaneous[1L], simultaneous_high = simultaneous[2L])
    }
    for (method in c("cusum", "support")) {
      yy <- dd[[paste0("reject_", method)]] - dd$reject_oracle
      n10 <- sum(yy == 1L); n01 <- sum(yy == -1L)
      component_error <- .025/(2*length(cells)*2)
      c10 <- cp_interval(n10,nr,component_error)
      c01 <- cp_interval(n01,nr,component_error)
      p10 <- cp_interval(n10,nr,.025)
      p01 <- cp_interval(n01,nr,.025)
      diff_rows[[length(diff_rows)+1L]] <- data.frame(
        cell = cell, contrast = paste0(method, "_minus_oracle"), trials = nr,
        discordant_10 = n10, discordant_01 = n01,
        estimate = mean(yy), MCSE = sd(yy)/sqrt(nr),
        pointwise_low = p10[1L]-p01[2L], pointwise_high = p10[2L]-p01[1L],
        simultaneous_low = c10[1L]-c01[2L],
        simultaneous_high = c10[2L]-c01[1L])
    }
    for (metric in c("date_error_cusum_over_k", "scaled_local_error", "date_error_support_over_k", "contamination_score_cusum",
                     "cusum_L2_residual", "cusum_sup_residual", "cusum_predicted_amplitude",
                     "support_L2_residual", "S_cusum_minus_oracle", "S_support_minus_oracle")) {
      yy <- dd[[metric]]
      se <- sd(yy)/sqrt(nr)
      half <- qt(.975, nr-1L)*se
      secondary_rows[[length(secondary_rows)+1L]] <- data.frame(
        cell = cell, metric = metric, estimate = mean(yy), SD = sd(yy), MCSE = se,
        exploratory_t_low = mean(yy)-half, exploratory_t_high = mean(yy)+half)
    }
    flags <- list(cusum_delayed = dd$cusum_date > dd$b, support_clipped = dd$support_clipped,
                  support_fallback = dd$support_fallback,
                  cusum_at_trim_boundary = dd$cusum_at_trim_boundary,
                  zero_joint_oracle = dd$joint_oracle == 0,
                  zero_joint_cusum = dd$joint_cusum == 0,
                  zero_joint_support = dd$joint_support == 0,
                  calibration_tie_oracle = dd$calibration_ties_oracle > 0,
                  calibration_tie_cusum = dd$calibration_ties_cusum > 0,
                  calibration_tie_support = dd$calibration_ties_support > 0)
    for (metric in names(flags)) {
      yy <- as.integer(flags[[metric]])
      ci <- cp_interval(sum(yy), nr, .05)
      frequency_rows[[length(frequency_rows)+1L]] <- data.frame(
        cell = cell, metric = metric, count = sum(yy), trials = nr, estimate = mean(yy),
        MCSE = sqrt(mean(yy)*(1-mean(yy))/nr), exploratory_CP95_low = ci[1L],
        exploratory_CP95_high = ci[2L])
    }
  }
  rejection <- do.call(rbind, rej_rows)
  paired <- do.call(rbind, diff_rows)
  secondary <- do.call(rbind, secondary_rows)
  write.csv(rejection, file.path(outdir, "rejection_summary.csv"), row.names = FALSE)
  write.csv(paired, file.path(outdir, "paired_summary.csv"), row.names = FALSE)
  write.csv(secondary, file.path(outdir, "exploratory_diagnostics.csv"), row.names = FALSE)
  write.csv(do.call(rbind,frequency_rows), file.path(outdir, "operational_frequencies.csv"),
            row.names = FALSE)
  plot_rejections(rejection, outdir)
  print(rejection, row.names = FALSE)
  print(paired, row.names = FALSE)
  cat("SYNTHETIC diagnostic; oracle critical values; no novelty or equivalence conclusion.\n")
  cat("All negative outcomes retained. Pointwise MCSE can be zero at boundary; inspect exact CI.\n")
  invisible(list(rejection = rejection, paired = paired, secondary = secondary))
}

plot_rejections <- function(tab,outdir) {
  pdf(file.path(outdir,"rejection_rates.pdf"),width=max(7.4, .85*length(unique(tab$cell))),height=4.8,useDingbats=FALSE)
  par(mar=c(5.2,4.2,1.6,1.0),las=1)
  cells <- unique(tab$cell)
  first <- tab[match(cells,tab$cell),]
  upper <- max(.15,tab$CI95_high)*1.08
  plot(NA,xlim=c(.5,length(cells)+.5),ylim=c(0,upper),xaxt="n",xlab="",
       ylab="Null rejection probability")
  axis(1,at=seq_along(cells),labels=paste0("n=",first$n,"\ng=",first$gamma),cex.axis=.8)
  abline(h=.05,lty=2,col="grey45")
  methods <- c("oracle","cusum","support")
  colors <- c("black","#B23A32","#245F99")
  for (m in seq_along(methods)) {
    z <- tab[tab$method==methods[m],];z <- z[match(cells,z$cell),]
    xx <- seq_along(cells)+(m-2)*.13
    arrows(xx,z$CI95_low,xx,z$CI95_high,angle=90,code=3,length=.035,col=colors[m])
    points(xx,z$estimate,pch=c(16,17,15)[m],col=colors[m])
  }
  legend("topright",legend=c("Known date","Raw-mean CUSUM","Known-endpoint date"),
         pch=c(16,17,15),col=colors,bty="n",cex=.8)
  mtext("Synthetic data; oracle Monte Carlo calibration; pointwise 95% exact intervals.",
        side=1,line=3.8,cex=.72)
  dev.off()
}


self_check <- function() {
  # Meaningful deterministic checks: independent exact integration and LS objective.
  for (n in 3:7) for (mask in 0:(2^n-1L)) {
    ind <- as.integer(intToBits(mask))[seq_len(n)]
    a <- sum(ind)
    if (a == 0) ref <- 0 else {
      ji <- c(0, cumsum(ind))[seq_len(n)]
      ref <- sum(vapply(seq_len(n), function(ii) {
        integrate(function(s) (ji[ii]-s*a)^2/a, (ii-1)/n, ii/n,
                  rel.tol = 1e-12)$value
      }, numeric(1L)))
    }
    stopifnot(abs(integral_statistic(ind)-ref) < 1e-10)
    dm <- drift_metrics(ind, ind, 1L, 1L, max(1L,floor(sqrt(n))), .1)
    stopifnot(dm[["L2_residual"]] == 0, dm[["sup_residual"]] == 0)
  }
  # This fixture separates unweighted CUSUM (j=4) from weighted LS (j=1).
  x <- c(2,0,0,0,-1,-1,-1,-1,-1,-1)
  jj <- trim_indices(length(x),.1)
  direct <- vapply(jj,function(j) j*(length(x)-j)/length(x)*
    (mean(x[seq_len(j)])-mean(x[seq.int(j+1L,length(x))])),numeric(1L))
  weighted <- direct^2/(jj*(1-jj/length(x)))
  stopifnot(cusum_date(x,.1)==4L,jj[which.max(abs(direct))]==4L,
            jj[which.max(weighted)]==1L,cusum_date(2*x+7,.1)==4L,
            cusum_date(rep(0,10),.1)==1L)
  # No random data in selfcheck. Exact location preserves within-segment ranks.
  uv <- cbind(u=c(.1,.9,.3,.8,.4,.7), v=c(.2,.8,.4,.6,.1,.9))
  loc <- location_from_uniforms(uv,.1,.5)
  stopifnot(identical(loc$x[1:3],uv[1:3,1]),
            max(abs(loc$x[4:6]-(uv[4:6,1]-.1/sqrt(2))))<1e-14,
            identical(loc$y,uv[,2]))
  base <- joint_indicators(uv[,1],uv[,2],3L,2L)
  stopifnot(identical(joint_indicators(loc$x,loc$y,3L,2L),base))
  stopifnot(identical(base, c(1L,0L,0L,0L,1L,0L))) # Actual cutoff inclusion.
  transformed <- joint_indicators(c(uv[1:3,1]^3, 5+exp(uv[4:6,1])),
                                    exp(uv[,2]),3L,2L)
  stopifnot(identical(base, transformed))
  # Y-only inclusion at i=2 must fail X, X-only at i=5 must fail Y.
  separated <- joint_indicators(uv[,1],c(.1,.2,.6,.9,.7,.8),3L,2L)
  stopifnot(identical(separated,c(1L,0L,0L,0L,0L,0L)))
  # An attainable exact segment cutoff: rank 1/(3+1) equals k/n=3/12.
  xedge <- c(.1,.3,.2,seq(1,9))
  yedge <- seq_len(12)
  stopifnot(1/(3+1) == 3/12,
            identical(joint_indicators(xedge,yedge,3L,3L),c(1L,rep(0L,11))))
  # Exact drift-integral vs independent quadrature for a nonzero signed difference.
  ie <- c(1L,1L,0L,0L,0L,1L); io <- c(0L,0L,0L,1L,0L,0L)
  dm <- drift_metrics(ie,io,4L,3L,2L,.2)
  cdf <- c(0,cumsum(ie-io)); n <- length(ie)
  amp <- .25*.2/sqrt(2)
  ref <- sum(vapply(seq_len(n), function(ii) integrate(function(s) {
    ((cdf[ii]-s*cdf[n+1L])/sqrt(2)+amp*(pmin(s/.5,1)-s))^2
  },(ii-1)/n,ii/n,rel.tol=1e-12)$value,numeric(1L)))
  stopifnot(cdf[n+1L] != 0,
            abs(dm[["L2_residual"]]^2-ref) < 1e-10,
            abs(dm[["sup_residual"]]-27/(20*sqrt(2))) < 1e-10)
  # Discordant-binomial component intervals agree with the independent R exact test.
  for (v in c(0L,1L,50L,100L)) {
    got <- cp_interval(v,100,.025)
    ref <- binom.test(v,100,conf.level=.975)$conf.int
    stopifnot(max(abs(got-ref))<1e-10)
  }
  cat("Deterministic integration, rank cutoffs, raw unweighted CUSUM, CI and drift checks passed.\n")
}


initialize_design <- function(profile,outdir) {
 cfg <- get_design(profile)
 stopifnot(!dir.exists(outdir));dir.create(outdir,recursive=TRUE)
 dput(cfg,file.path(outdir,"executed_design.R"))
 writeLines(capture.output(sessionInfo()),file.path(outdir,"session_info.txt"))
 writeLines(c(R.version.string,paste("R.home:",R.home()),
  paste("commandArgs:",paste(commandArgs(),collapse=" "))),file.path(outdir,"runtime.txt"))
 cells <- expand.grid(n=cfg$n,gamma=cfg$gamma,KEEP.OUT.ATTRS=FALSE)
 RNGkind("L'Ecuyer-CMRG","Inversion","Rejection");set.seed(cfg$seed)
 total <- nrow(cells)*cfg$R
 cal <- eva <- matrix(0L,nrow=7L,ncol=total)
 stream <- .Random.seed
 for(id in seq_len(total)) {
  cal[,id] <- stream;eva[,id] <- parallel::nextRNGStream(stream)
  stream <- parallel::nextRNGStream(eva[,id])
 }
 saveRDS(list(profile=profile,design=cfg,calibration=cal,evaluation=eva),file.path(outdir,"seeds.rds"),version=3)
 records <- data.frame(cell=rep(rep(seq_len(nrow(cells)),each=cfg$R),each=2L),
  replication=rep(rep(seq_len(cfg$R),times=nrow(cells)),each=2L),
  role=rep(c("calibration","evaluation"),times=total))
 states <- matrix(0L,nrow=2L*total,ncol=7L)
 states[seq.int(1L,2L*total,by=2L),] <- t(cal)
 states[seq.int(2L,2L*total,by=2L),] <- t(eva)
 colnames(states) <- paste0("seed",seq_len(7L))
 write.csv(cbind(records,states),file.path(outdir,"rng_streams.csv"),row.names=FALSE)
 writeLines(capture.output(RNGkind()),file.path(outdir,"rng_kind.txt"))
 cat("Initialized",profile,"with",total,"evaluation stream assignments.\n")
}

run_chunk <- function(profile,cell,start,count,seedfile,outdir) {
 cfg <- get_design(profile)
 cells <- expand.grid(n=cfg$n,gamma=cfg$gamma,KEEP.OUT.ATTRS=FALSE)
 stopifnot(cell>=1L,cell<=nrow(cells),start>=1L,count>=1L,count<=10L,
           start+count-1L<=cfg$R,!dir.exists(outdir))
 seeds <- readRDS(seedfile)
 stopifnot(identical(seeds$profile,profile),identical(seeds$design,cfg),
           identical(dim(seeds$calibration),c(7L,as.integer(nrow(cells)*cfg$R))),
           identical(dim(seeds$evaluation),dim(seeds$calibration)))
 dir.create(outdir,recursive=TRUE)
 n <- as.integer(cells$n[cell]);gamma <- cells$gamma[cell]
 k <- as.integer(floor(sqrt(n)));b <- as.integer(floor(cfg$tau*n))
 RNGkind("L'Ecuyer-CMRG","Inversion","Rejection")
 for(replication in seq.int(start,start+count-1L)) {
  id <- (cell-1L)*cfg$R+replication
  assign(".Random.seed",seeds$calibration[,id],envir=.GlobalEnv)
  calibration <- numeric(cfg$B)
  for(bb in seq_len(cfg$B)) {
   uv <- clayton_one(n)
   calibration[bb] <- integral_statistic(joint_indicators(uv[,1L],uv[,2L],b,k))
   write_append(data.frame(cell=cell,replication=replication,calibration_id=bb,
    statistic=calibration[bb]),file.path(outdir,"calibration_values.csv"))
  }
  assign(".Random.seed",seeds$evaluation[,id],envir=.GlobalEnv)
  sam <- make_sample(n,gamma,cfg$tau)
  jb <- cusum_date(sam$x,cfg$trim);js <- support_date(sam$x,cfg$trim)
  io <- joint_indicators(sam$x,sam$y,b,k)
  ib <- joint_indicators(sam$x,sam$y,jb,k)
  it <- joint_indicators(sam$x,sam$y,js$j,k)
  values <- c(oracle=integral_statistic(io),cusum=integral_statistic(ib),support=integral_statistic(it))
  pv <- vapply(values,function(x)(1+sum(calibration>=x))/(cfg$B+1),numeric(1L))
  mb <- drift_metrics(ib,io,jb,b,k,sam$Delta);mt <- drift_metrics(it,io,js$j,b,k,sam$Delta)
  row <- data.frame(cell=cell,n=n,gamma=gamma,replication=replication,k=k,Delta=sam$Delta,
   b=b,cusum_date=jb,support_date=js$j,support_clipped=js$clipped,support_fallback=js$fallback,
   cusum_at_trim_boundary=jb%in%range(trim_indices(n,cfg$trim)),
   date_error_cusum_over_k=(jb-b)/k,scaled_local_error=3*gamma^2*(jb-b)/k,
   date_error_support_over_k=(js$j-b)/k,contamination_score_cusum=max(jb-b,0)*sam$Delta/sqrt(k),
   cusum_L2_residual=unname(mb["L2_residual"]),cusum_sup_residual=unname(mb["sup_residual"]),
   cusum_predicted_amplitude=unname(mb["predicted_amplitude"]),support_L2_residual=unname(mt["L2_residual"]),
   S_oracle=values[["oracle"]],S_cusum=values[["cusum"]],S_support=values[["support"]],
   S_cusum_minus_oracle=values[["cusum"]]-values[["oracle"]],S_support_minus_oracle=values[["support"]]-values[["oracle"]],
   p_oracle=pv[["oracle"]],p_cusum=pv[["cusum"]],p_support=pv[["support"]],
   reject_oracle=as.integer(pv[["oracle"]]<=cfg$alpha),reject_cusum=as.integer(pv[["cusum"]]<=cfg$alpha),
   reject_support=as.integer(pv[["support"]]<=cfg$alpha),
   joint_oracle=sum(io),joint_cusum=sum(ib),joint_support=sum(it),
   calibration_ties_oracle=sum(calibration==values[["oracle"]]),
   calibration_ties_cusum=sum(calibration==values[["cusum"]]),
   calibration_ties_support=sum(calibration==values[["support"]]))
  write_append(row,file.path(outdir,"trials.csv"))
  cat("cell",cell,"replication",replication,"completed\n");flush.console()
 }
 writeLines("complete",file.path(outdir,"COMPLETE.txt"))
}

aggregate_design <- function(profile,listfile,outdir) {
 cfg <- get_design(profile)
 stopifnot(!dir.exists(outdir));dir.create(outdir,recursive=TRUE)
 files <- readLines(listfile,warn=FALSE)
 dat <- do.call(rbind,lapply(files,read.csv,check.names=FALSE))
 cells <- expand.grid(n=cfg$n,gamma=cfg$gamma,KEEP.OUT.ATTRS=FALSE)
 expected <- expand.grid(replication=seq_len(cfg$R),cell=seq_len(nrow(cells)))
 key <- function(x)paste(x$cell,x$replication,sep=":")
 stopifnot(nrow(dat)==nrow(expected),!anyDuplicated(key(dat)),setequal(key(dat),key(expected)),
  all(dat$n==cells$n[dat$cell]),all(dat$gamma==cells$gamma[dat$cell]))
 dat <- dat[order(dat$cell,dat$replication),]
 write.csv(dat,file.path(outdir,"trials.csv"),row.names=FALSE)
 saveRDS(list(design=cfg,data_label="SYNTHETIC",trials=dat),file.path(outdir,"complete_results.rds"),version=3)
 dput(cfg,file.path(outdir,"executed_design.R"))
 summarize_results(dat,cfg,outdir)
}

main <- function() {
 a <- commandArgs(trailingOnly=TRUE)
 receipt <- Sys.getenv("ACADEMIC_RUN_RECEIPT");token <- Sys.getenv("ACADEMIC_RUN_TOKEN")
 if(!nzchar(token)||!file.exists(receipt)||!identical(readLines(receipt,n=1L,warn=FALSE),token))
  stop("Use the bounded run_full.py launcher; it records seeds, sources and process-tree time limits.")
 pid <- Sys.getenv("ACADEMIC_R_PID_FILE");if(nzchar(pid))writeLines(as.character(Sys.getpid()),pid)
 if(length(a)==1L&&a[1L]=="selfcheck")self_check()
 else if(length(a)==3L&&a[1L]=="initialize")initialize_design(a[2L],a[3L])
 else if(length(a)==7L&&a[1L]=="chunk")run_chunk(a[2L],as.integer(a[3L]),as.integer(a[4L]),
  as.integer(a[5L]),a[6L],a[7L])
 else if(length(a)==4L&&a[1L]=="aggregate")aggregate_design(a[2L],a[3L],a[4L])
 else stop("Use the bounded launcher: run_full.py smoke|full [--resume RUN_FOLDER]")
}
if(sys.nframe()==0L)main()
