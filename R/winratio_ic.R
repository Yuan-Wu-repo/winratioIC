# winratio_ic.R
# User-facing companion implementation for the proposed win-ratio methodology.
#
# Supported observation schemes in this first version:
#   1) margin1 = interval, margin2 = interval
#   2) margin1 = right,    margin2 = interval
#
# Here margin1 is always the higher-priority outcome.
# "right" means an exactly observed event time subject to ordinary right censoring.
#
# The internal interval/interval engine is extracted from exponential+BB1(7).txt.
# The internal right/interval engine is extracted from exponential+BB1(9).txt.
# Simulation-specific data-generation and Monte Carlo reporting code are intentionally omitted.
#
# The public tau_l argument controls the lower integration bound. For
# interval/interval, both prioritized components are integrated over [tau_l, tau].
# For right/interval, the exactly observed/right-censored primary component remains
# integrated over [0, tau], while the interval-censored secondary component uses
# [tau_l, tau]. tau_l = 0 computes the untrimmed exploratory target.
#
# IMPORTANT for the mixed terminal-outcome setting:
# For margin2, the supplied interval must already represent the actually observed
# secondary-outcome information after any truncation at the end of primary follow-up.
# The analysis function does not reconstruct a visit history or re-apply simulation recoding.
# Optional visit_lower/visit_upper columns may be supplied for an interval margin
# when exact reproduction of the original simulation knot placement is desired.

`%||%` <- function(x, y) if (is.null(x)) y else x

.wr_require_columns <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing)) {
    stop("Missing required data column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.wr_validate_margin <- function(margin, which_margin) {
  if (!is.list(margin) || is.null(margin$type)) {
    stop(which_margin, " must be a list with a 'type' element.", call. = FALSE)
  }
  type <- tolower(as.character(margin$type)[1])
  if (!type %in% c("interval", "right")) {
    stop(which_margin, "$type must be 'interval' or 'right'.", call. = FALSE)
  }
  if (type == "interval") {
    if (is.null(margin$left) || is.null(margin$right)) {
      stop(which_margin, " with type='interval' must specify 'left' and 'right' column names.", call. = FALSE)
    }
  } else {
    if (is.null(margin$time) || is.null(margin$status)) {
      stop(which_margin, " with type='right' must specify 'time' and 'status' column names.", call. = FALSE)
    }
  }
  type
}

.wr_prepare_interval <- function(data, margin, label) {
  .wr_require_columns(data, c(margin$left, margin$right))
  L <- as.numeric(data[[margin$left]])
  R <- as.numeric(data[[margin$right]])

  if (length(L) != nrow(data) || length(R) != nrow(data)) {
    stop(label, ": interval columns have invalid lengths.", call. = FALSE)
  }
  if (anyNA(L) || anyNA(R)) {
    stop(label, ": missing interval endpoints are not allowed.", call. = FALSE)
  }
  if (any(!is.finite(L)) || any(L < 0)) {
    stop(label, ": left endpoints must be finite and >= 0.", call. = FALSE)
  }
  if (any(R <= 0, na.rm = TRUE)) {
    stop(label, ": finite right endpoints must be > 0.", call. = FALSE)
  }
  if (any(is.finite(R) & L > R)) {
    stop(label, ": each interval must satisfy left <= right.", call. = FALSE)
  }
  if (any(is.finite(R) & L == R)) {
    stop(label, ": exact observations (left == right) are not supported for type='interval'.", call. = FALSE)
  }

  # Public representation is the usual (L, R] form:
  #   left censored:  L=0, R<Inf
  #   interval:       0<L<R<Inf
  #   right censored: L>=0, R=Inf
  left_cens  <- (L == 0 & is.finite(R))
  int_cens   <- (L > 0 & is.finite(R))
  right_cens <- is.infinite(R) & R > 0

  if (any(!(left_cens | int_cens | right_cens))) {
    stop(label, ": intervals must be left-, interval-, or right-censored in (L,R] notation.", call. = FALSE)
  }

  # Convert to the simulation-code representation.
  #
  # By default, the ordinary observed censoring interval (L,R] supplies all
  # likelihood-relevant endpoints. For left- and right-censored observations,
  # the non-likelihood endpoint is duplicated so that spline evaluation remains
  # finite. This is the natural user-facing mode when only censoring intervals
  # are available.
  #
  # OPTIONAL exact-reproduction mode:
  # If both margin$visit_lower and margin$visit_upper are supplied, those
  # columns are used as the finite inspection pair (u,v) passed to the original
  # estimator. The censoring category is STILL determined only by (L,R].
  # This allows exact reproduction of the submitted simulation implementation,
  # whose data-adaptive spline-knot rule used both finite inspection times even
  # when only one of them entered a left/right-censoring likelihood term.
  has_vlo <- !is.null(margin$visit_lower)
  has_vup <- !is.null(margin$visit_upper)
  if (xor(has_vlo, has_vup)) {
    stop(label, ": supply both 'visit_lower' and 'visit_upper', or neither.", call. = FALSE)
  }

  if (has_vlo && has_vup) {
    .wr_require_columns(data, c(margin$visit_lower, margin$visit_upper))
    tu <- as.numeric(data[[margin$visit_lower]])
    tv <- as.numeric(data[[margin$visit_upper]])

    if (anyNA(tu) || anyNA(tv) || any(!is.finite(tu)) || any(!is.finite(tv)) ||
        any(tu < 0) || any(tv < 0) || any(tu > tv)) {
      stop(label, ": visit_lower/visit_upper must be finite, nonnegative, and satisfy visit_lower <= visit_upper.", call. = FALSE)
    }

    tol <- 1e-8 * pmax(1, abs(L), ifelse(is.finite(R), abs(R), 0), abs(tu), abs(tv))
    inconsistent <-
      (left_cens  & abs(tu - R) > tol) |
      (int_cens   & (abs(tu - L) > tol | abs(tv - R) > tol)) |
      (right_cens & abs(tv - L) > tol)

    if (any(inconsistent)) {
      stop(label, ": visit_lower/visit_upper are inconsistent with the supplied (L,R] censoring interval for at least one row.", call. = FALSE)
    }
  } else {
    tu <- numeric(length(L))
    tv <- numeric(length(L))

    tu[left_cens] <- R[left_cens]
    tv[left_cens] <- R[left_cens]

    tu[int_cens] <- L[int_cens]
    tv[int_cens] <- R[int_cens]

    tu[right_cens] <- L[right_cens]
    tv[right_cens] <- L[right_cens]
  }

  list(
    tu = tu,
    tv = tv,
    delta1 = as.numeric(left_cens),
    delta2 = as.numeric(int_cens),
    delta3 = as.numeric(right_cens)
  )
}

.wr_prepare_right <- function(data, margin, label) {
  .wr_require_columns(data, c(margin$time, margin$status))
  time <- as.numeric(data[[margin$time]])
  status <- data[[margin$status]]

  if (anyNA(time) || any(!is.finite(time)) || any(time < 0)) {
    stop(label, ": time must be finite and >= 0.", call. = FALSE)
  }

  if (is.logical(status)) status <- as.integer(status)
  status <- suppressWarnings(as.numeric(status))
  if (anyNA(status) || any(!status %in% c(0, 1))) {
    stop(label, ": status must be coded 1=event, 0=right censored.", call. = FALSE)
  }

  list(
    time = time,
    delta_event = status,
    delta_censor = 1 - status
  )
}

.wr_split_treatment <- function(data, treatment) {
  if (!is.character(treatment) || length(treatment) != 1L) {
    stop("'treatment' must be the name of one column in data.", call. = FALSE)
  }
  .wr_require_columns(data, treatment)
  z <- data[[treatment]]
  if (anyNA(z)) stop("Treatment variable contains missing values.", call. = FALSE)

  if (is.logical(z)) z <- as.integer(z)

  if (is.numeric(z) || is.integer(z)) {
    u <- sort(unique(as.numeric(z)))
    if (!identical(u, c(0, 1))) {
      stop("Numeric treatment must be coded 0=control and 1=treatment.", call. = FALSE)
    }
    return(list(control = which(z == 0), treatment = which(z == 1), control_label = 0, treatment_label = 1))
  }

  if (is.factor(z)) {
    lev <- levels(z)
    lev <- lev[lev %in% unique(as.character(z))]
    if (length(lev) != 2L) stop("Factor treatment must have exactly two observed levels.", call. = FALSE)
    return(list(
      control = which(as.character(z) == lev[1]),
      treatment = which(as.character(z) == lev[2]),
      control_label = lev[1], treatment_label = lev[2]
    ))
  }

  # Recognize common character labels, otherwise require an ordered factor or 0/1 coding.
  zz <- tolower(trimws(as.character(z)))
  u <- unique(zz)
  if (setequal(u, c("control", "treatment"))) {
    return(list(
      control = which(zz == "control"), treatment = which(zz == "treatment"),
      control_label = unique(as.character(z)[zz == "control"])[1],
      treatment_label = unique(as.character(z)[zz == "treatment"])[1]
    ))
  }
  stop("Character treatment labels are ambiguous. Recode to 0/1 or use a factor whose first level is control and second level is treatment.", call. = FALSE)
}

.wr_subset_prepared <- function(x, idx) lapply(x, function(v) v[idx])


.wr_ii_engine <- local({
  size <- 1L
  tau_analysis <- 5
  k.i <- 6L
  sk.i <- k.i^2

  gg_pre <- function(qqu_1,qqu_2,qqv_1,qqv_2,
                     data1,data2,data3,data4,data5,data6,data7,data8,data9,
                     noze1,noze2,noze3,noze4,noze5,noze6,noze7,noze8,noze9,
                     sum11,sum22,sum33,sum44,sum55,sum66,sum77,sum88,sum99,
                     w)
  {
    ddpre <- matrix(0, (sk.i+2*k.i), 1)

    idx1 <- which(data1!=0); idx2 <- which(data2!=0); idx3 <- which(data3!=0)
    idx4 <- which(data4!=0); idx5 <- which(data5!=0); idx6 <- which(data6!=0)
    idx7 <- which(data7!=0); idx8 <- which(data8!=0); idx9 <- which(data9!=0)

    # Weighted score contributions: sum_i w_i * (1/p_i) * dp_i/dtheta
    mat1 <- qqu_1[,idx1,drop=FALSE] %*% (t(qqu_2[,idx1,drop=FALSE]) * (w[idx1] / sum11[idx1]))
    mat2 <- qqu_1[,idx2,drop=FALSE] %*% (t(qqv_2[,idx2,drop=FALSE]) * (w[idx2] / sum22[idx2])) -
      qqu_1[,idx2,drop=FALSE] %*% (t(qqu_2[,idx2,drop=FALSE]) * (w[idx2] / sum22[idx2]))
    mat3 <- qqu_1[,idx3,drop=FALSE] %*% (matrix(1,noze3,k.i) * (w[idx3] / sum33[idx3])) -
      qqu_1[,idx3,drop=FALSE] %*% (t(qqv_2[,idx3,drop=FALSE]) * (w[idx3] / sum33[idx3]))
    mat4 <- qqv_1[,idx4,drop=FALSE] %*% (t(qqu_2[,idx4,drop=FALSE]) * (w[idx4] / sum44[idx4])) -
      qqu_1[,idx4,drop=FALSE] %*% (t(qqu_2[,idx4,drop=FALSE]) * (w[idx4] / sum44[idx4]))
    mat5 <- qqv_1[,idx5,drop=FALSE] %*% (t(qqv_2[,idx5,drop=FALSE]) * (w[idx5] / sum55[idx5])) -
      qqv_1[,idx5,drop=FALSE] %*% (t(qqu_2[,idx5,drop=FALSE]) * (w[idx5] / sum55[idx5])) -
      qqu_1[,idx5,drop=FALSE] %*% (t(qqv_2[,idx5,drop=FALSE]) * (w[idx5] / sum55[idx5])) +
      qqu_1[,idx5,drop=FALSE] %*% (t(qqu_2[,idx5,drop=FALSE]) * (w[idx5] / sum55[idx5]))
    mat6 <- qqv_1[,idx6,drop=FALSE] %*% (matrix(1,noze6,k.i) * (w[idx6] / sum66[idx6])) -
      qqv_1[,idx6,drop=FALSE] %*% (t(qqv_2[,idx6,drop=FALSE]) * (w[idx6] / sum66[idx6])) -
      qqu_1[,idx6,drop=FALSE] %*% (matrix(1,noze6,k.i) * (w[idx6] / sum66[idx6])) +
      qqu_1[,idx6,drop=FALSE] %*% (t(qqv_2[,idx6,drop=FALSE]) * (w[idx6] / sum66[idx6]))
    mat7 <- (matrix(1,k.i,noze7) - qqv_1[,idx7,drop=FALSE]) %*% (t(qqu_2[,idx7,drop=FALSE]) * (w[idx7] / sum77[idx7]))
    mat8 <- matrix(1,k.i,noze8) %*% (t(qqv_2[,idx8,drop=FALSE]) * (w[idx8] / sum88[idx8])) -
      matrix(1,k.i,noze8) %*% (t(qqu_2[,idx8,drop=FALSE]) * (w[idx8] / sum88[idx8])) -
      qqv_1[,idx8,drop=FALSE] %*% (t(qqv_2[,idx8,drop=FALSE]) * (w[idx8] / sum88[idx8])) +
      qqv_1[,idx8,drop=FALSE] %*% (t(qqu_2[,idx8,drop=FALSE]) * (w[idx8] / sum88[idx8]))
    mat9 <- -matrix(1,k.i,noze9) %*% (t(qqv_2[,idx9,drop=FALSE]) * (w[idx9] / sum99[idx9])) -
      qqv_1[,idx9,drop=FALSE] %*% (matrix(1,noze9,k.i) * (w[idx9] / sum99[idx9])) +
      qqv_1[,idx9,drop=FALSE] %*% (t(qqv_2[,idx9,drop=FALSE]) * (w[idx9] / sum99[idx9]))

    for (jj in 1:sk.i)
    {
      if(as.integer(jj/k.i)==jj/k.i) { m <- jj/k.i; n <- k.i }
      else { m <- as.integer(jj/k.i)+1; n <- jj-as.integer(jj/k.i)*k.i }
      ddpre[jj,1] <- mat1[m,n]+mat2[m,n]+mat3[m,n]+mat4[m,n]+mat5[m,n]+mat6[m,n]+mat7[m,n]+mat8[m,n]+mat9[m,n]
    }

    vec7 <- matrix(1,1,noze7) %*% (t(qqu_2[,idx7,drop=FALSE]) * (w[idx7] / sum77[idx7]))
    vec8 <- matrix(1,1,noze8) %*% (t(qqv_2[,idx8,drop=FALSE]) * (w[idx8] / sum88[idx8])) -
      matrix(1,1,noze8) %*% (t(qqu_2[,idx8,drop=FALSE]) * (w[idx8] / sum88[idx8]))
    vec9 <- -matrix(1,1,noze9) %*% (t(qqv_2[,idx9,drop=FALSE]) * (w[idx9] / sum99[idx9]))

    for (jj in (sk.i+1):(sk.i+k.i))
    {
      n <- jj-sk.i
      ddpre[jj,1] <- vec7[n]+vec8[n]+vec9[n]
    }

    vec3 <- qqu_1[,idx3,drop=FALSE] %*% (matrix(1,noze3,1) * (w[idx3] / sum33[idx3]))
    vec6 <- qqv_1[,idx6,drop=FALSE] %*% (matrix(1,noze6,1) * (w[idx6] / sum66[idx6])) -
      qqu_1[,idx6,drop=FALSE] %*% (matrix(1,noze6,1) * (w[idx6] / sum66[idx6]))
    vec9 <- -qqv_1[,idx9,drop=FALSE] %*% (matrix(1,noze9,1) * (w[idx9] / sum99[idx9]))

    for (jj in (sk.i+k.i+1):(sk.i+2*k.i))
    {
      m <- jj-(sk.i+k.i)
      ddpre[jj,1] <- vec3[m]+vec6[m]+vec9[m]
    }

    ddpre
  }


  hh_pre <- function(qqu_1,qqu_2,qqv_1,qqv_2,
                     data1,data2,data3,data4,data5,data6,data7,data8,data9,
                     sum11,sum22,sum33,sum44,sum55,sum66,sum77,sum88,sum99,
                     mat11,mat22,mat33,mat44,mat55,mat66,mat77,mat88,mat99,
                     vec77,vec88,vec99,vec33,vec66,vec99new,
                     w)
  {
    hhpre <- matrix(0, (sk.i+2*k.i), (sk.i+2*k.i))

    idx1 <- which(data1!=0); idx2 <- which(data2!=0); idx3 <- which(data3!=0)
    idx4 <- which(data4!=0); idx5 <- which(data5!=0); idx6 <- which(data6!=0)
    idx7 <- which(data7!=0); idx8 <- which(data8!=0); idx9 <- which(data9!=0)

    # Weighted Hessian: - sum_i w_i * (grad p_i)(grad p_i)^T / p_i^2
    ww11 <- w[idx1]/(sum11[idx1]^2)
    ww22 <- w[idx2]/(sum22[idx2]^2)
    ww33 <- w[idx3]/(sum33[idx3]^2)
    ww44 <- w[idx4]/(sum44[idx4]^2)
    ww55 <- w[idx5]/(sum55[idx5]^2)
    ww66 <- w[idx6]/(sum66[idx6]^2)
    ww77 <- w[idx7]/(sum77[idx7]^2)
    ww88 <- w[idx8]/(sum88[idx8]^2)
    ww99 <- w[idx9]/(sum99[idx9]^2)

    mat111 <- mat11 %*% (t(mat11) * ww11)
    mat222 <- mat22 %*% (t(mat22) * ww22)
    mat333 <- mat33 %*% (t(mat33) * ww33)
    mat444 <- mat44 %*% (t(mat44) * ww44)
    mat555 <- mat55 %*% (t(mat55) * ww55)
    mat666 <- mat66 %*% (t(mat66) * ww66)
    mat777 <- mat77 %*% (t(mat77) * ww77)
    mat888 <- mat88 %*% (t(mat88) * ww88)
    mat999 <- mat99 %*% (t(mat99) * ww99)

    for(ii in 1:sk.i)
      for(jj in ii:sk.i)
        hhpre[ii,jj] <- -mat111[ii,jj]-mat222[ii,jj]-mat333[ii,jj]-mat444[ii,jj]-mat555[ii,jj]-mat666[ii,jj]-mat777[ii,jj]-mat888[ii,jj]-mat999[ii,jj]

    vec777 <- mat77 %*% (t(vec77) * ww77)
    vec888 <- mat88 %*% (t(vec88) * ww88)
    vec999 <- mat99 %*% (t(vec99) * ww99)

    for (ii in 1:sk.i)
      for (jj in (sk.i+1):(sk.i+k.i))
        hhpre[ii,jj] <- -vec777[ii, jj-sk.i]-vec888[ii, jj-sk.i]-vec999[ii, jj-sk.i]

    vec333 <- mat33 %*% (t(vec33) * ww33)
    vec666 <- mat66 %*% (t(vec66) * ww66)
    vec999b <- mat99 %*% (t(vec99new) * ww99)

    for (ii in 1:sk.i)
      for (jj in (sk.i+k.i+1):(sk.i+2*k.i))
        hhpre[ii,jj] <- -vec333[ii, jj-sk.i-k.i]-vec666[ii, jj-sk.i-k.i]-vec999b[ii, jj-sk.i-k.i]

    vec7772 <- vec77 %*% (t(vec77) * ww77)
    vec8882 <- vec88 %*% (t(vec88) * ww88)
    vec9992 <- vec99 %*% (t(vec99) * ww99)

    for (ii in (sk.i+1):(sk.i+k.i))
      for (jj in ii:(sk.i+k.i))
        hhpre[ii,jj] <- -vec7772[ii-sk.i, jj-sk.i]-vec8882[ii-sk.i, jj-sk.i]-vec9992[ii-sk.i, jj-sk.i]

    vec9993 <- vec99 %*% (t(vec99new) * ww99)
    for (ii in (sk.i+1):(sk.i+k.i))
      for (jj in (sk.i+k.i+1):(sk.i+2*k.i))
        hhpre[ii,jj] <- -vec9993[ii-sk.i, jj-sk.i-k.i]

    vec3333 <- vec33 %*% (t(vec33) * ww33)
    vec6663 <- vec66 %*% (t(vec66) * ww66)
    vec9994 <- vec99new %*% (t(vec99new) * ww99)

    for (ii in (sk.i+k.i+1):(sk.i+2*k.i))
      for (jj in ii:(sk.i+2*k.i))
        hhpre[ii,jj] <- -vec3333[ii-sk.i-k.i, jj-sk.i-k.i]-vec6663[ii-sk.i-k.i, jj-sk.i-k.i]-vec9994[ii-sk.i-k.i, jj-sk.i-k.i]

    # symmetrize
    hhpre <- hhpre + t(hhpre)
    for (ii in 1:(sk.i+2*k.i)) hhpre[ii,ii] <- hhpre[ii,ii]/2
    hhpre
  }


  get_knots <- function(tu_1,tu_2,tv_1,tv_2,
                        delta1_1,delta2_1,delta3_1,
                        delta1_2,delta2_2,delta3_2)
  {
    err <- 0.5

    knot_pre_1 <- c(tu_1[delta1_1!=0], tu_1[delta2_1!=0],
                    tv_1[delta2_1!=0], tv_1[delta3_1!=0])
    knot_pre_2 <- c(tu_2[delta1_2!=0], tu_2[delta2_2!=0],
                    tv_2[delta2_2!=0], tv_2[delta3_2!=0])

    ## Original rule has k.i - 3 interior knots at fractions
    ## 1/(k.i-2), ..., (k.i-3)/(k.i-2).
    probs <- seq_len(k.i - 3L) / (k.i - 2L)

    interior1 <- as.numeric(quantile(knot_pre_1, probs = probs, names = FALSE))
    interior2 <- as.numeric(quantile(knot_pre_2, probs = probs, names = FALSE))

    ## FIX: boundary must contain every x passed to splineDesign.
    upper1 <- max(c(tu_1, tv_1, tau_analysis), na.rm = TRUE) + err
    upper2 <- max(c(tu_2, tv_2, tau_analysis), na.rm = TRUE) + err

    knot1 <- c(rep(0, 3), interior1, rep(upper1, 3))
    knot2 <- c(rep(0, 3), interior2, rep(upper2, 3))

    stopifnot(length(knot1) == k.i + 3L,
              length(knot2) == k.i + 3L)

    list(knot1 = knot1, knot2 = knot2)
  }


  sumvectorfunction <- function(alpha,qqu_1,qqu_2,qqv_1,qqv_2)
  {
    sumuu=rep(0,size)
    for(i in 1:k.i) for(j in 1:k.i) sumuu=sumuu+alpha[(i-1)*k.i+j,1]*qqu_1[i,]*qqu_2[j,]
    sumuv=rep(0,size)
    for(i in 1:k.i) for(j in 1:k.i) sumuv=sumuv+alpha[(i-1)*k.i+j,1]*qqu_1[i,]*qqv_2[j,]
    sumvu=rep(0,size)
    for(i in 1:k.i) for(j in 1:k.i) sumvu=sumvu+alpha[(i-1)*k.i+j,1]*qqv_1[i,]*qqu_2[j,]
    sumvv=rep(0,size)
    for(i in 1:k.i) for(j in 1:k.i) sumvv=sumvv+alpha[(i-1)*k.i+j,1]*qqv_1[i,]*qqv_2[j,]

    sumu2=rep(0,size)
    for(i in 1:k.i) for(j in 1:k.i) sumu2=sumu2+alpha[(i-1)*k.i+j,1]*qqu_2[j,]
    sumv2=rep(0,size)
    for(i in 1:k.i) for(j in 1:k.i) sumv2=sumv2+alpha[(i-1)*k.i+j,1]*qqv_2[j,]
    sumu1=rep(0,size)
    for(i in 1:k.i) for(j in 1:k.i) sumu1=sumu1+alpha[(i-1)*k.i+j,1]*qqu_1[i,]
    sumv1=rep(0,size)
    for(i in 1:k.i) for(j in 1:k.i) sumv1=sumv1+alpha[(i-1)*k.i+j,1]*qqv_1[i,]

    sumu2m=rep(0,size)
    for(i in 1:k.i) sumu2m=sumu2m+alpha[(sk.i)+i,1]*qqu_2[i,]
    sumv2m=rep(0,size)
    for(i in 1:k.i) sumv2m=sumv2m+alpha[(sk.i)+i,1]*qqv_2[i,]
    sumu1m=rep(0,size)
    for(i in 1:k.i) sumu1m=sumu1m+alpha[(sk.i+k.i)+i,1]*qqu_1[i,]
    sumv1m=rep(0,size)
    for(i in 1:k.i) sumv1m=sumv1m+alpha[(sk.i+k.i)+i,1]*qqv_1[i,]

    sum11<-sumuu
    sum22<-sumuv-sumuu
    sum33<-(sumu1+sumu1m)-sumuv
    sum44<-sumvu-sumuu
    sum55<-sumvv-sumvu-sumuv+sumuu
    sum66<-(sumv1+sumv1m)-sumvv-(sumu1+sumu1m)+sumuv
    sum77<-(sumu2+sumu2m)-sumvu
    sum88<-(sumv2+sumv2m)-(sumu2+sumu2m)-sumvv+sumvu
    sum99<-1-(sumv2+sumv2m)-(sumv1+sumv1m)+sumvv

    list(sum11=sum11,sum22=sum22,sum33=sum33,sum44=sum44,sum55=sum55,sum66=sum66,sum77=sum77,sum88=sum88,sum99=sum99)
  }


  sumfunction <- function(sum11,sum22,sum33,sum44,sum55,sum66,sum77,sum88,sum99,
                          data1,data2,data3,data4,data5,data6,data7,data8,data9,
                          w)
  {
    probs <- list(sum11,sum22,sum33,sum44,sum55,sum66,sum77,sum88,sum99)
    dats  <- list(data1,data2,data3,data4,data5,data6,data7,data8,data9)

    total <- 0
    for (j in seq_along(probs)) {
      idx <- dats[[j]] != 0
      if (!any(idx)) next
      pj <- probs[[j]][idx]
      if (any(!is.finite(pj)) || any(pj <= 0)) return(-Inf)
      total <- total + sum(w[idx] * log(pj))
    }
    total
  }


  totest_diag <- function(tu_1,tu_2,tv_1,tv_2,
                     delta1_1,delta2_1,delta3_1,delta1_2,delta2_2,delta3_2,
                     w = rep(1, size),
                     alpha_init = NULL,
                     max_halvings = 10L)
  {
    knot12 <- get_knots(tu_1,tu_2,tv_1,tv_2,delta1_1,delta2_1,delta3_1,delta1_2,delta2_2,delta3_2)
    knot1 <- knot12$knot1
    knot2 <- knot12$knot2

    # I-splines via B-splines
    bsp <- t(splines::splineDesign(knots = c(0, knot1, max(knot1)), x = tu_1, ord = 4))
    pre.isp <- apply(bsp, 2, rev); isp <- apply(pre.isp, 2, cumsum); isp <- apply(isp, 2 ,rev); isp <- isp[-1,]
    qqu_1 <- isp

    bsp <- t(splines::splineDesign(knots = c(0, knot1, max(knot1)), x = tv_1, ord = 4))
    pre.isp <- apply(bsp, 2, rev); isp <- apply(pre.isp, 2, cumsum); isp <- apply(isp, 2 ,rev); isp <- isp[-1,]
    qqv_1 <- isp

    bsp <- t(splines::splineDesign(knots = c(0, knot2, max(knot2)), x = tu_2, ord = 4))
    pre.isp <- apply(bsp, 2, rev); isp <- apply(pre.isp, 2, cumsum); isp <- apply(isp, 2 ,rev); isp <- isp[-1,]
    qqu_2 <- isp

    bsp <- t(splines::splineDesign(knots = c(0, knot2, max(knot2)), x = tv_2, ord = 4))
    pre.isp <- apply(bsp, 2, rev); isp <- apply(pre.isp, 2, cumsum); isp <- apply(isp, 2 ,rev); isp <- isp[-1,]
    qqv_2 <- isp

    data1 <- delta1_1*delta1_2
    data2 <- delta1_1*delta2_2
    data3 <- delta1_1*delta3_2
    data4 <- delta2_1*delta1_2
    data5 <- delta2_1*delta2_2
    data6 <- delta2_1*delta3_2
    data7 <- delta3_1*delta1_2
    data8 <- delta3_1*delta2_2
    data9 <- delta3_1*delta3_2

    noze1 <- sum(data1); noze2 <- sum(data2); noze3 <- sum(data3)
    noze4 <- sum(data4); noze5 <- sum(data5); noze6 <- sum(data6)
    noze7 <- sum(data7); noze8 <- sum(data8); noze9 <- sum(data9)

    mat11 <- matrix(0,sk.i,noze1); mat22 <- matrix(0,sk.i,noze2); mat33 <- matrix(0,sk.i,noze3)
    mat44 <- matrix(0,sk.i,noze4); mat55 <- matrix(0,sk.i,noze5); mat66 <- matrix(0,sk.i,noze6)
    mat77 <- matrix(0,sk.i,noze7); mat88 <- matrix(0,sk.i,noze8); mat99 <- matrix(0,sk.i,noze9)

    for (ii in 1:sk.i)
    {
      if(as.integer(ii/k.i)==ii/k.i) { m1=ii/k.i; n1=k.i }
      else { m1=as.integer(ii/k.i)+1; n1=ii-as.integer(ii/k.i)*k.i }

      mat11[ii,] <- qqu_1[m1,data1!=0]*qqu_2[n1,data1!=0]
      mat22[ii,] <- qqu_1[m1,data2!=0]*qqv_2[n1,data2!=0]-qqu_1[m1,data2!=0]*qqu_2[n1,data2!=0]
      mat33[ii,] <- qqu_1[m1,data3!=0]-qqu_1[m1,data3!=0]*qqv_2[n1,data3!=0]
      mat44[ii,] <- qqv_1[m1,data4!=0]*qqu_2[n1,data4!=0]-qqu_1[m1,data4!=0]*qqu_2[n1,data4!=0]
      mat55[ii,] <- qqv_1[m1,data5!=0]*qqv_2[n1,data5!=0]-qqv_1[m1,data5!=0]*qqu_2[n1,data5!=0]-
        qqu_1[m1,data5!=0]*qqv_2[n1,data5!=0]+qqu_1[m1,data5!=0]*qqu_2[n1,data5!=0]
      mat66[ii,] <- qqv_1[m1,data6!=0]-qqv_1[m1,data6!=0]*qqv_2[n1,data6!=0]-qqu_1[m1,data6!=0]+qqu_1[m1,data6!=0]*qqv_2[n1,data6!=0]
      mat77[ii,] <- qqu_2[n1,data7!=0]-qqv_1[m1,data7!=0]*qqu_2[n1,data7!=0]
      mat88[ii,] <- qqv_2[n1,data8!=0]-qqu_2[n1,data8!=0]-qqv_1[m1,data8!=0]*qqv_2[n1,data8!=0]+qqv_1[m1,data8!=0]*qqu_2[n1,data8!=0]
      mat99[ii,] <- -qqv_2[n1,data9!=0]-qqv_1[m1,data9!=0]+qqv_1[m1,data9!=0]*qqv_2[n1,data9!=0]
    }

    vec77 <- qqu_2[,data7!=0,drop=FALSE]
    vec88 <- qqv_2[,data8!=0,drop=FALSE]-qqu_2[,data8!=0,drop=FALSE]
    vec99 <- -qqv_2[,data9!=0,drop=FALSE]

    vec33 <- qqu_1[,data3!=0,drop=FALSE]
    vec66 <- qqv_1[,data6!=0,drop=FALSE]-qqu_1[,data6!=0,drop=FALSE]
    vec99new <- -qqv_1[,data9!=0,drop=FALSE]

    ## Starting value.  The original fit uses the original generic interior
    ## start.  Weighted-bootstrap refits may instead be warm-started from a
    ## slightly interiorized version of the original unweighted MLE.
    p_alpha <- sk.i + 2*k.i
    if (is.null(alpha_init)) {
      alpha <- matrix(0.5/p_alpha, p_alpha, 1)
    } else {
      avec <- as.numeric(alpha_init)
      if (length(avec) != p_alpha || any(!is.finite(avec)))
        stop("Invalid alpha_init")
      alpha <- matrix(avec, p_alpha, 1)
      ## Numerical safeguard only: keep the supplied start feasible.
      alpha[alpha < 0] <- 0
      if (sum(alpha) >= 1) alpha <- alpha / sum(alpha) * (1 - 1e-10)
    }

    sumvector <- sumvectorfunction(alpha,qqu_1,qqu_2,qqv_1,qqv_2)
    sum11 <- sumvector$sum11; sum22 <- sumvector$sum22; sum33 <- sumvector$sum33
    sum44 <- sumvector$sum44; sum55 <- sumvector$sum55; sum66 <- sumvector$sum66
    sum77 <- sumvector$sum77; sum88 <- sumvector$sum88; sum99 <- sumvector$sum99

    summl1 <- sumfunction(sum11,sum22,sum33,sum44,sum55,sum66,sum77,sum88,sum99,
                          data1,data2,data3,data4,data5,data6,data7,data8,data9,
                          w)

    dd <- gg_pre(qqu_1,qqu_2,qqv_1,qqv_2,
                 data1,data2,data3,data4,data5,data6,data7,data8,data9,
                 noze1,noze2,noze3,noze4,noze5,noze6,noze7,noze8,noze9,
                 sum11,sum22,sum33,sum44,sum55,sum66,sum77,sum88,sum99,
                 w)

    hhh <- hh_pre(qqu_1,qqu_2,qqv_1,qqv_2,
                  data1,data2,data3,data4,data5,data6,data7,data8,data9,
                  sum11,sum22,sum33,sum44,sum55,sum66,sum77,sum88,sum99,
                  mat11,mat22,mat33,mat44,mat55,mat66,mat77,mat88,mat99,
                  vec77,vec88,vec99,vec33,vec66,vec99new,
                  w)

    # original generalized gradient-projection algorithm with implementation safeguards
    index=0
    rowa1<--1
    rowa=0
    whilecheck<-0
    kk=-1
    kk1=1
    kk2=1
    stop=0
    fail_reason <- "success"
    fail_detail <- NA_character_

    while(rowa1!=rowa)
    {
      if (whilecheck==0)
      {
        whilecheck<-1
        hreg <- -hhh + diag(0.01, (sk.i+2*k.i))
        kk1 <- rcond(hreg)
        if (!is.finite(kk1) || kk1 < 10^(-15)) {
          fail_reason <- "hessian_ill_conditioned"
          fail_detail <- paste0("initial rcond=", kk1)
          stop <- 1
          break
        }
        ww <- solve(hreg)
        gt=ww%*%dd
        rowa=0
        pw=diag(1,(sk.i+2*k.i))
        ddd=pw%*%gt
        for(i in 1:(sk.i+2*k.i)) if (abs(ddd[i])<10^(-7)/5) ddd[i]=0
      } else {
        ## FIX 2a: aa exists only when the active set is nonempty.
        if (rowa > 0) {
          amat <- aa %*% ww %*% t(aa)
          kk2 <- rcond(amat)
          if (!is.finite(kk2) || kk2 < 10^(-15)) {
            fail_reason <- "active_set_ill_conditioned"
            fail_detail <- paste0("rcond=", kk2)
            stop <- 1
            break
          }
          pw <- diag(1,(sk.i+2*k.i)) -
                ww %*% t(aa) %*% solve(amat) %*% aa
        } else {
          pw <- diag(1,(sk.i+2*k.i))
        }
        ddd=pw%*%gt
        for(i in 1:(sk.i+2*k.i)) if (abs(ddd[i])<10^(-7)/5) ddd[i]=0
      }

      while(max(abs(ddd))>0.0001)
      {
        iiii=0; tep=0; tep1=0
        sumc=sum(alpha); sumd=sum(ddd)

        if(1-sumc>0.0000001&sumd>0.0000001)
        {
          iiii=iiii+1
          tep[iiii]=(1-sumc)/sumd
          tep1[iiii]=-1
        }

        for( i in 1:(k.i))
        {
          for( j in 1:(k.i))
          {
            if (alpha[(i-1)*k.i+j]>0.0000001&ddd[(i-1)*k.i+j]<0)
            {
              iiii=iiii+1
              tep[iiii]=(0.0000001-alpha[(i-1)*k.i+j])/ddd[(i-1)*k.i+j]
              tep1[iiii]=(i-1)*k.i+j
            }
          }
        }

        for(i in (sk.i+1):(sk.i+2*k.i))
        {
          if(alpha[i]>0&ddd[i]<0)
          {
            iiii=iiii+1
            tep[iiii]=-alpha[i]/ddd[i]
            tep1[iiii]=i
          }
        }

        rrtt=tep[1]; kkkk=1
        if(iiii>0) { rrtt=min(tep); kkkk=which.min(tep) }
        rrtt=min(rrtt,0.5)

        summl2=summl1-1; kk=-1
        while(summl2<=summl1 & kk<=max_halvings)
        {
          kk=kk+1
          alphat=alpha+0.5^kk*rrtt*ddd
          if (sum(alphat)>1+10^(-5)) { fail_reason <- "mass_constraint_violation"; fail_detail <- paste0("sum(alphat)=", sum(alphat)); stop=1; break }

          u1=tep1[kkkk]; k5=0; add=0
          if(u1!=-1&kk==0&tep[kkkk]<=0.5)
          {
            if(u1>sk.i) alphat[u1]=0 else alphat[u1]=0.0000001

            for(i in 1:(sk.i+2*k.i))
            {
              if(((alphat[i]<10^(-8))||(abs(alphat[i]-10^(-7))<10^(-8)&i<=sk.i))&i!=tep1[kkkk])
              {
                if(!any(index==i)) { k5=k5+1; add[k5]=i }
              }
            }
          }

          for(i in 1:(sk.i+2*k.i)) if(alphat[i]<10^(-8)) alphat[i]=0

          sumvector <- sumvectorfunction(alphat,qqu_1,qqu_2,qqv_1,qqv_2)
          sum11 <- sumvector$sum11; sum22 <- sumvector$sum22; sum33 <- sumvector$sum33
          sum44 <- sumvector$sum44; sum55 <- sumvector$sum55; sum66 <- sumvector$sum66
          sum77 <- sumvector$sum77; sum88 <- sumvector$sum88; sum99 <- sumvector$sum99

          summl2 <- sumfunction(sum11,sum22,sum33,sum44,sum55,sum66,sum77,sum88,sum99,
                                data1,data2,data3,data4,data5,data6,data7,data8,data9,
                                w)

          if(summl2==Inf) { fail_reason <- "trial_loglik_posInf"; fail_detail <- "summl2 == Inf"; stop=1; break }
        }
        if (stop==1) break
        if (kk>max_halvings) {
          fail_reason <- "line_search_failure"
          fail_detail <- paste0("No likelihood increase after >", max_halvings, " step halvings")
          stop=1
          break
        }

        if(kk==0&tep[kkkk]<=0.5)
        {
          rowa=rowa+1
          index[rowa]=tep1[kkkk]
          anew=matrix(0,1,(sk.i+2*k.i))

          if(tep1[kkkk]==(-1)) { for(i in 1:(sk.i+2*k.i)) anew[1,i]=1 }
          else { anew[1,tep1[kkkk]]=-1 }

          if(rowa==1) aa=anew
          if(rowa>1) aa=rbind(aa,anew)

          if(k5!=0)
          {
            for(i in 1:k5) index[rowa+i]=add[i]
            rowa=rowa+k5

            anew=matrix(0,k5,(sk.i+2*k.i))
            for(i in 1:k5) anew[i,add[i]]=-1
            aa=rbind(aa,anew)
          }
        }

        alpha=alphat

        sumvector <- sumvectorfunction(alpha,qqu_1,qqu_2,qqv_1,qqv_2)
        sum11 <- sumvector$sum11; sum22 <- sumvector$sum22; sum33 <- sumvector$sum33
        sum44 <- sumvector$sum44; sum55 <- sumvector$sum55; sum66 <- sumvector$sum66
        sum77 <- sumvector$sum77; sum88 <- sumvector$sum88; sum99 <- sumvector$sum99

        summl1 <- sumfunction(sum11,sum22,sum33,sum44,sum55,sum66,sum77,sum88,sum99,
                              data1,data2,data3,data4,data5,data6,data7,data8,data9,
                              w)

        dd <- gg_pre(qqu_1,qqu_2,qqv_1,qqv_2,
                     data1,data2,data3,data4,data5,data6,data7,data8,data9,
                     noze1,noze2,noze3,noze4,noze5,noze6,noze7,noze8,noze9,
                     sum11,sum22,sum33,sum44,sum55,sum66,sum77,sum88,sum99,
                     w)

        hhh <- hh_pre(qqu_1,qqu_2,qqv_1,qqv_2,
                      data1,data2,data3,data4,data5,data6,data7,data8,data9,
                      sum11,sum22,sum33,sum44,sum55,sum66,sum77,sum88,sum99,
                      mat11,mat22,mat33,mat44,mat55,mat66,mat77,mat88,mat99,
                      vec77,vec88,vec99,vec33,vec66,vec99new,
                      w)

        kk1 <- rcond(-hhh+diag(0.01, (sk.i+2*k.i)))
        if(kk1<10^(-15)) { fail_reason <- "hessian_ill_conditioned"; fail_detail <- paste0("rcond=", kk1); stop=1; break }
        ww <- solve(-hhh+diag(0.01, (sk.i+2*k.i)))

        gt <- ww%*%dd

        ## FIX 2b: do not evaluate an active-set matrix when rowa == 0.
        if (rowa > 0) {
          amat <- aa %*% ww %*% t(aa)
          kk2 <- rcond(amat)
          if (!is.finite(kk2) || kk2 < 10^(-15)) {
            fail_reason <- "active_set_ill_conditioned"
            fail_detail <- paste0("rcond=", kk2)
            stop <- 1
            break
          }
          pw <- diag(1,(sk.i+2*k.i)) -
                ww %*% t(aa) %*% solve(amat) %*% aa
        } else {
          pw <- diag(1,(sk.i+2*k.i))
        }
        ddd <- pw%*%gt
        for(i in 1:(sk.i+2*k.i)) if (abs(ddd[i])<10^(-7)/5) ddd[i]=0
      }

      if (stop==1) break

      ## FIX 2c: Lagrange multipliers are defined only for a nonempty
      ## active set.  If rowa == 0 there is no aa and no multiplier to compute.
      rowa1 <- rowa
      if (rowa > 0) {
        amat <- aa %*% ww %*% t(aa)
        kk2 <- rcond(amat)
        if (!is.finite(kk2) || kk2 < 10^(-15)) {
          fail_reason <- "active_set_ill_conditioned"
          fail_detail <- paste0("rcond=", kk2)
          stop <- 1
          break
        }
        lambda <- solve(amat) %*% aa %*% gt

        if (rowa > 1 && min(lambda) < 0)
        {
          rowa=rowa-1
          ik=which.min(lambda)
          index=index[-ik]
          if (rowa>1) aa=aa[-ik,]
          else { if(rowa==1) aa=t(aa[-ik,]) }
        }
      }
    }

    if (stop==0) {
      list(alpha = alpha, knot1 = knot1, knot2 = knot2,
           status = "success", detail = NA_character_)
    } else {
      list(alpha = -9999, knot1 = -9999, knot2 = -9999,
           status = fail_reason, detail = fail_detail)
    }
  }


  mar.surv1 <- function(t1, knot11, alpha.joint.1, alpha.mar.1.1)
  {
    knotb <- c(0, knot11, max(knot11))
    bsp <- t(splines::splineDesign(knots = knotb, x = t1, ord = 4))
    pre.isp <- rev(bsp); isp <- cumsum(pre.isp); isp <- rev(isp); isp <- isp[-1]
    q1 <- isp

    tep_matrix <- matrix(c(alpha.joint.1),k.i,k.i)
    tep_q1 <- matrix(q1, 1, k.i)
    sum1 <- tep_q1 %*% t(tep_matrix)
    sum2 <- alpha.mar.1.1 * tep_q1
    1 - (sum(sum1) + sum(sum2))
  }


  mar.dcdf0 <- function(t1, knot01, alpha.joint.0, alpha.mar.0.1)
  {
    knotb <- knot01
    ord <- 3
    bsp <- t(splines::splineDesign(knots = knotb, x = t1, ord = ord))
    msp <- bsp / (knotb[(1 + ord) : (k.i + ord)] - knotb[1 : k.i]) * 3
    q1 <- msp

    tep_matrix <- matrix(c(alpha.joint.0),k.i,k.i)
    tep_q1 <- matrix(q1, 1, k.i)
    sum1 <- tep_q1 %*% t(tep_matrix)
    sum2 <- alpha.mar.0.1 * tep_q1
    sum(sum1) + sum(sum2)
  }


  mar.multip <- function(t, knot01, knot11, alpha.joint.0, alpha.mar.0.1, alpha.joint.1, alpha.mar.1.1)
  {
    mar.surv1(t, knot11, alpha.joint.1, alpha.mar.1.1) * mar.dcdf0(t, knot01, alpha.joint.0, alpha.mar.0.1)
  }


  est.win.ratio1 <- function(knot01, knot11, alpha.joint.0, alpha.mar.0.1,
                             alpha.joint.1, alpha.mar.1.1, tau_l = 0)
  {
    ## Restricted primary comparison over [tau_l, tau_analysis].
    start <- tau_l
    end <- min(max(knot01), max(knot11), tau_analysis)
    if (end <= start) stop("Invalid primary integral boundaries")
    parts <- 1000L
    dif <- (end - start) / parts
    integrate.vec <- vapply(seq_len(parts), function(ind) {
      mar.multip(start + ind*dif,
                 knot01, knot11, alpha.joint.0, alpha.mar.0.1,
                 alpha.joint.1, alpha.mar.1.1) * dif
    }, numeric(1))
    sum(integrate.vec)
  }


  ## Evaluate the I-spline CDF basis at a scalar time, consistently with the likelihood.
  get_I_basis <- function(t, knot) {
    bsp <- as.numeric(splines::splineDesign(
      knots = c(0, knot, max(knot)), x = t, ord = 4
    ))
    rev(cumsum(rev(bsp)))[-1L]
  }


  ## Evaluate the corresponding M-spline density basis.
  get_M_basis <- function(t, knot) {
    ord <- 3L
    bsp <- as.numeric(splines::splineDesign(knots = knot, x = t, ord = ord))
    3 * bsp / (knot[(1L + ord):(k.i + ord)] - knot[1L:k.i])
  }


  ## Evaluate primary survival at the finite analysis horizon tau_analysis.
  joint.surv1 <- function(t2, knot11, knot12, alpha.joint.1,
                          alpha.mar.1.1, alpha.mar.1.2) {
    q1 <- get_I_basis(tau_analysis, knot11)
    q2 <- get_I_basis(t2, knot12)
    A <- matrix(as.numeric(alpha.joint.1), k.i, k.i)
    F1 <- sum((colSums(A) + alpha.mar.1.1) * q1)
    F2 <- sum((rowSums(A) + alpha.mar.1.2) * q2)
    F12 <- sum(A * outer(q2, q1))
    1 - F1 - F2 + F12
  }


  ## partial_2 S_0(tau_analysis,t2), which is nonpositive.
  joint.dcdf0 <- function(t2, knot01, knot02, alpha.joint.0,
                          alpha.mar.0.2) {
    q1 <- get_I_basis(tau_analysis, knot01)
    m2 <- get_M_basis(t2, knot02)
    A <- matrix(as.numeric(alpha.joint.0), k.i, k.i)
    f2 <- sum((rowSums(A) + alpha.mar.0.2) * m2)
    dF12 <- sum(A * outer(m2, q1))
    dF12 - f2
  }


  joint.multip <- function(t, knot01, knot02, knot11, knot12,
                           alpha.joint.0, alpha.mar.0.2,
                           alpha.joint.1, alpha.mar.1.1, alpha.mar.1.2) {
    -joint.surv1(t, knot11, knot12, alpha.joint.1,
                 alpha.mar.1.1, alpha.mar.1.2) *
      joint.dcdf0(t, knot01, knot02, alpha.joint.0, alpha.mar.0.2)
  }


  est.win.ratio2 <- function(knot01, knot02, knot11, knot12,
                             alpha.joint.0, alpha.mar.0.2,
                             alpha.joint.1, alpha.mar.1.1, alpha.mar.1.2,
                             tau_l = 0) {
    ## Restricted secondary comparison over [tau_l, tau_analysis].
    start <- tau_l
    end <- min(max(knot02), max(knot12), tau_analysis)
    if (end <= start) stop("Invalid secondary integral boundaries")
    parts <- 1000L
    dif <- (end - start) / parts
    integrate.vec <- vapply(seq_len(parts), function(ind) {
      joint.multip(start + ind*dif, knot01, knot02, knot11, knot12,
                   alpha.joint.0, alpha.mar.0.2,
                   alpha.joint.1, alpha.mar.1.1, alpha.mar.1.2) * dif
    }, numeric(1))
    sum(integrate.vec)
  }


  est_win_prob <- function(spl0, spl1, tau_l = 0) {
    knot01 <- spl0$knot1; knot02 <- spl0$knot2
    knot11 <- spl1$knot1; knot12 <- spl1$knot2

    alpha.joint.0 <- spl0$alpha[1:(k.i*k.i)]
    alpha.joint.1 <- spl1$alpha[1:(k.i*k.i)]

    alpha.mar.0.2 <- spl0$alpha[(k.i*k.i + 1):(k.i*k.i + k.i)]
    alpha.mar.0.1 <- spl0$alpha[(k.i*k.i + k.i + 1):(k.i*k.i + 2*k.i)]

    alpha.mar.1.2 <- spl1$alpha[(k.i*k.i + 1):(k.i*k.i + k.i)]
    alpha.mar.1.1 <- spl1$alpha[(k.i*k.i + k.i + 1):(k.i*k.i + 2*k.i)]

    est.win.ratio1(knot01, knot11, alpha.joint.0, alpha.mar.0.1,
                   alpha.joint.1, alpha.mar.1.1, tau_l = tau_l) +
      est.win.ratio2(knot01, knot02, knot11, knot12,
                     alpha.joint.0, alpha.mar.0.2,
                     alpha.joint.1, alpha.mar.1.1, alpha.mar.1.2,
                     tau_l = tau_l)
  }


  est_log_wr <- function(spl0, spl1, eps = 1e-12, tau_l = 0) {
    win  <- est_win_prob(spl0, spl1, tau_l = tau_l)
    lose <- est_win_prob(spl1, spl0, tau_l = tau_l)
    if (!is.finite(win) || !is.finite(lose) || win <= eps || lose <= eps) return(NA_real_)
    log(win/lose)
  }

  validate_fit <- function(ans, tol = 1e-7) {
    if (is.null(ans$status)) {
      ans$status <- "invalid_fit_object"
      ans$detail <- "Missing status"
      return(ans)
    }

    if (!identical(ans$status, "success")) return(ans)

    p <- sk.i + 2L*k.i
    a <- as.numeric(ans$alpha)

    if (length(a) != p || any(!is.finite(a))) {
      ans$status <- "invalid_coefficient_vector"
      ans$detail <- "Coefficient vector has wrong length or non-finite entries"
      return(ans)
    }

    if (any(a < -tol)) {
      ans$status <- "nonnegativity_constraint_failure"
      ans$detail <- paste0("minimum alpha=", min(a))
      return(ans)
    }

    if (sum(a) > 1 + tol) {
      ans$status <- "mass_constraint_failure"
      ans$detail <- paste0("sum(alpha)=", sum(a))
      return(ans)
    }

    if (any(!is.finite(ans$knot1)) || any(!is.finite(ans$knot2))) {
      ans$status <- "invalid_knot_vector"
      ans$detail <- "Non-finite knot value"
      return(ans)
    }

    ans
  }


  safe_fit <- function(dat, w = rep(1, size),
                       alpha_init = NULL,
                       max_halvings = 10L) {

    warning_messages <- character(0)

    ans <- tryCatch(
      withCallingHandlers(
        totest_diag(
          dat$tu_1, dat$tu_2, dat$tv_1, dat$tv_2,
          dat$delta1_1, dat$delta2_1, dat$delta3_1,
          dat$delta1_2, dat$delta2_2, dat$delta3_2,
          w = w,
          alpha_init = alpha_init,
          max_halvings = max_halvings
        ),
        warning = function(wrn) {
          warning_messages <<- c(warning_messages, conditionMessage(wrn))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        list(
          alpha = -9999,
          knot1 = -9999,
          knot2 = -9999,
          status = "R_error",
          detail = conditionMessage(e)
        )
      }
    )

    ans <- validate_fit(ans)

    ans$n_warnings <- length(warning_messages)
    ans$warning_detail <- if (length(warning_messages) == 0L) {
      NA_character_
    } else {
      paste(unique(warning_messages), collapse = " | ")
    }

    ans
  }


  est_win_components_diag <- function(spl0, spl1, tau_l = 0) {

    knot01 <- spl0$knot1
    knot02 <- spl0$knot2
    knot11 <- spl1$knot1
    knot12 <- spl1$knot2

    alpha.joint.0 <- spl0$alpha[1:(k.i*k.i)]
    alpha.joint.1 <- spl1$alpha[1:(k.i*k.i)]

    alpha.mar.0.2 <-
      spl0$alpha[(k.i*k.i + 1):(k.i*k.i + k.i)]
    alpha.mar.0.1 <-
      spl0$alpha[(k.i*k.i + k.i + 1):(k.i*k.i + 2*k.i)]

    alpha.mar.1.2 <-
      spl1$alpha[(k.i*k.i + 1):(k.i*k.i + k.i)]
    alpha.mar.1.1 <-
      spl1$alpha[(k.i*k.i + k.i + 1):(k.i*k.i + 2*k.i)]

    primary <- tryCatch(
      est.win.ratio1(
        knot01, knot11,
        alpha.joint.0, alpha.mar.0.1,
        alpha.joint.1, alpha.mar.1.1,
        tau_l = tau_l
      ),
      error = function(e) e
    )

    if (inherits(primary, "error")) {
      return(list(
        status = "primary_functional_error",
        detail = conditionMessage(primary),
        primary = NA_real_,
        secondary = NA_real_,
        total = NA_real_
      ))
    }

    if (length(primary) != 1L || !is.finite(primary)) {
      return(list(
        status = "primary_functional_nonfinite",
        detail = NA_character_,
        primary = as.numeric(primary)[1],
        secondary = NA_real_,
        total = NA_real_
      ))
    }

    secondary <- tryCatch(
      est.win.ratio2(
        knot01, knot02, knot11, knot12,
        alpha.joint.0, alpha.mar.0.2,
        alpha.joint.1, alpha.mar.1.1, alpha.mar.1.2,
        tau_l = tau_l
      ),
      error = function(e) e
    )

    if (inherits(secondary, "error")) {
      return(list(
        status = "secondary_functional_error",
        detail = conditionMessage(secondary),
        primary = primary,
        secondary = NA_real_,
        total = NA_real_
      ))
    }

    if (length(secondary) != 1L || !is.finite(secondary)) {
      return(list(
        status = "secondary_functional_nonfinite",
        detail = NA_character_,
        primary = primary,
        secondary = as.numeric(secondary)[1],
        total = NA_real_
      ))
    }

    total <- primary + secondary

    if (!is.finite(total)) {
      return(list(
        status = "total_winprob_nonfinite",
        detail = NA_character_,
        primary = primary,
        secondary = secondary,
        total = total
      ))
    }

    list(
      status = "success",
      detail = NA_character_,
      primary = primary,
      secondary = secondary,
      total = total
    )
  }


  est_log_wr_diag <- function(spl0, spl1, eps = 1e-12, tau_l = 0) {

    win <- est_win_components_diag(spl0, spl1, tau_l = tau_l)

    if (win$status != "success") {
      return(list(
        status = paste0("win_", win$status),
        detail = win$detail,
        logwr = NA_real_,
        win = NA_real_,
        lose = NA_real_
      ))
    }

    lose <- est_win_components_diag(spl1, spl0, tau_l = tau_l)

    if (lose$status != "success") {
      return(list(
        status = paste0("loss_", lose$status),
        detail = lose$detail,
        logwr = NA_real_,
        win = win$total,
        lose = NA_real_
      ))
    }

    if (!is.finite(win$total) || win$total <= eps) {
      return(list(
        status = "nonpositive_or_nonfinite_win",
        detail = NA_character_,
        logwr = NA_real_,
        win = win$total,
        lose = lose$total
      ))
    }

    if (!is.finite(lose$total) || lose$total <= eps) {
      return(list(
        status = "nonpositive_or_nonfinite_loss",
        detail = NA_character_,
        logwr = NA_real_,
        win = win$total,
        lose = lose$total
      ))
    }

    lw <- log(win$total / lose$total)

    if (!is.finite(lw)) {
      return(list(
        status = "nonfinite_logWR",
        detail = NA_character_,
        logwr = lw,
        win = win$total,
        lose = lose$total
      ))
    }

    list(
      status = "success",
      detail = NA_character_,
      logwr = lw,
      win = win$total,
      lose = lose$total
    )
  }


  make_interior_start <- function(alpha_hat, mix = 1e-4) {

    p <- sk.i + 2L*k.i
    generic <- matrix(0.5/p, p, 1)

    avec <- as.numeric(alpha_hat)

    if (length(avec) != p || any(!is.finite(avec))) {
      stop("Invalid original alpha_hat")
    }

    ah <- matrix(avec, p, 1)
    ah[ah < 0] <- 0

    if (sum(ah) >= 1) {
      ah <- ah / sum(ah) * (1 - 1e-10)
    }

    (1 - mix) * ah + mix * generic
  }


  fit_bootstrap_with_rescue <- function(dat, w, alpha_hat) {

    ## First attempt: start very close to the unweighted MLE.
    start1 <- make_interior_start(alpha_hat, mix = 1e-4)

    fit1 <- safe_fit(
      dat,
      w = w,
      alpha_init = start1,
      max_halvings = 30L
    )

    if (identical(fit1$status, "success")) {
      return(list(
        fit = fit1,
        initial_status = fit1$status,
        initial_detail = fit1$detail,
        rescue_used = FALSE,
        rescue_status = NA_character_,
        rescue_detail = NA_character_
      ))
    }

    ## Rescue attempt: SAME weights, more interior start, longer line search.
    start2 <- make_interior_start(alpha_hat, mix = 0.05)

    fit2 <- safe_fit(
      dat,
      w = w,
      alpha_init = start2,
      max_halvings = 60L
    )

    list(
      fit = fit2,
      initial_status = fit1$status,
      initial_detail = fit1$detail,
      rescue_used = TRUE,
      rescue_status = fit2$status,
      rescue_detail = fit2$detail
    )
  }


  configure <- function(n, k, tau_value) {
    size <<- as.integer(n)
    k.i <<- as.integer(k)
    sk.i <<- k.i^2
    tau_analysis <<- as.numeric(tau_value)
    invisible(NULL)
  }

  fit_original <- function(dat, k, tau_value) {
    n <- length(dat$tu_1)
    configure(n, k, tau_value)
    fit <- safe_fit(
      dat,
      w = rep(1, n),
      alpha_init = NULL,
      max_halvings = 10L
    )
    list(fit = fit)
  }

  fit_bootstrap <- function(dat, w, alpha_hat, k, tau_value) {
    configure(length(dat$tu_1), k, tau_value)
    fit_bootstrap_with_rescue(dat, w, alpha_hat)
  }

  evaluate <- function(fit0, fit1, k, tau_value, tau_l = 0) {
    configure(1L, k, tau_value)
    z <- est_log_wr_diag(fit0, fit1, tau_l = tau_l)
    if (!identical(z$status, "success")) return(z)
    list(
      status = "success",
      detail = z$detail,
      logwr = z$logwr,
      win = z$win,
      lose = z$lose,
      wr = exp(z$logwr)
    )
  }

  list(
    fit_original = fit_original,
    fit_bootstrap = fit_bootstrap,
    evaluate = evaluate
  )
})


.wr_ri_engine <- local({
  size <- 1L
  tau <- 5
  k.i <- 6L
  n.spline <- k.i
  sk.i <- k.i^2
  spline.ord <- 4L
  p_alpha <- sk.i + 2L*k.i
  dt <- 0.02
  ds <- 0.02

  sumvectorfunction <- function(alpha, ispx, mspx, ispu, ispv)
  {
    mu <- t(matrix(alpha[1 : sk.i], k.i, k.i))
    omega <- matrix(alpha[(sk.i+k.i+1) : (sk.i+2*k.i)], k.i, 1)
    pi <- matrix(alpha[(sk.i+1) : (sk.i+k.i)], k.i, 1)

    fxu <- diag(t(ispx) %*% mu %*% ispu)
    dfxu <- diag(t(mspx) %*% mu %*% ispu)
    fxv <- diag(t(ispx) %*% mu %*% ispv)
    dfxv <- diag(t(mspx) %*% mu %*% ispv)

    fx <- diag(t(ispx) %*% mu %*% matrix(1, n.spline, size)) + c(t(ispx) %*% omega)
    dfx <- diag(t(mspx) %*% mu %*% matrix(1, n.spline, size)) + c(t(mspx) %*% omega)
    fu <- diag(matrix(1, size, n.spline) %*% mu %*% ispu) + c(t(ispu) %*% pi)
    fv <- diag(matrix(1, size, n.spline) %*% mu %*% ispv) + c(t(ispv) %*% pi)

    sum11 <- dfxu
    sum22 <- dfxv - dfxu
    sum33 <- dfx - dfxv
    sum44 <- fu - fxu
    sum55 <- fv - fu - fxv + fxu
    sum66 <- 1 - fv - fx + fxv

    list(sum11=sum11, sum22=sum22, sum33=sum33, sum44=sum44, sum55=sum55, sum66=sum66)
  }


  sumfunction_w <- function(sum11,sum22,sum33,sum44,sum55,sum66,
                            data1,data2,data3,data4,data5,data6,
                            w)
  {
    probs <- list(sum11,sum22,sum33,sum44,sum55,sum66)
    dats  <- list(data1,data2,data3,data4,data5,data6)

    total <- 0
    for (j in seq_along(probs)) {
      idx <- dats[[j]] != 0
      if (!any(idx)) next

      pj <- probs[[j]][idx]

      ## Numerical safeguard only: candidate points outside the valid
      ## likelihood region are rejected by the line search.
      if (any(!is.finite(pj)) || any(pj <= 0)) return(-Inf)

      total <- total + sum(w[idx] * log(pj))
    }

    total
  }


  gg_pre_w <- function(ispx, mspx, ispu, ispv,
                       data1, data2, data3, data4, data5, data6,
                       noze1, noze2, noze3, noze4, noze5, noze6,
                       sum11, sum22, sum33, sum44, sum55, sum66,
                       w)
  {
    ddpre <- matrix(0, (sk.i+2*k.i), 1)

    idx1 <- which(data1 != 0); idx2 <- which(data2 != 0); idx3 <- which(data3 != 0)
    idx4 <- which(data4 != 0); idx5 <- which(data5 != 0); idx6 <- which(data6 != 0)

    wrow <- function(A, ww) sweep(A, 1, ww, `*`)

    mat1 <- matrix(0, k.i, k.i)
    if (length(idx1)) {
      A <- t(ispu[, idx1, drop=FALSE]) / sum11[idx1]
      A <- wrow(A, w[idx1])
      mat1 <- mspx[, idx1, drop=FALSE] %*% A
    }

    mat2 <- matrix(0, k.i, k.i)
    if (length(idx2)) {
      A1 <- t(ispv[, idx2, drop=FALSE]) / sum22[idx2]
      A2 <- t(ispu[, idx2, drop=FALSE]) / sum22[idx2]
      A  <- wrow(A1 - A2, w[idx2])
      mat2 <- mspx[, idx2, drop=FALSE] %*% A
    }

    mat3 <- matrix(0, k.i, k.i)
    if (length(idx3)) {
      A1 <- matrix(1, length(idx3), k.i) / sum33[idx3]
      A2 <- t(ispv[, idx3, drop=FALSE]) / sum33[idx3]
      A  <- wrow(A1 - A2, w[idx3])
      mat3 <- mspx[, idx3, drop=FALSE] %*% A
    }

    mat4 <- matrix(0, k.i, k.i)
    if (length(idx4)) {
      B <- t(ispu[, idx4, drop=FALSE]) / sum44[idx4]
      B <- wrow(B, w[idx4])
      mat4 <- (matrix(1, k.i, length(idx4)) - ispx[, idx4, drop=FALSE]) %*% B
    }

    mat5 <- matrix(0, k.i, k.i)
    if (length(idx5)) {
      Bv <- t(ispv[, idx5, drop=FALSE]) / sum55[idx5]
      Bu <- t(ispu[, idx5, drop=FALSE]) / sum55[idx5]
      B  <- Bv - Bu
      B  <- wrow(B, w[idx5])

      term1 <- matrix(1, k.i, length(idx5)) %*% B
      term2 <- - ispx[, idx5, drop=FALSE] %*% wrow(Bv, w[idx5])
      term3 <- + ispx[, idx5, drop=FALSE] %*% wrow(Bu, w[idx5])
      mat5 <- term1 + term2 + term3
    }

    mat6 <- matrix(0, k.i, k.i)
    if (length(idx6)) {
      Bv <- t(ispv[, idx6, drop=FALSE]) / sum66[idx6]
      Bv <- wrow(Bv, w[idx6])

      term1 <- -matrix(1, k.i, length(idx6)) %*% Bv
      term2 <- -ispx[, idx6, drop=FALSE] %*% wrow(matrix(1, length(idx6), k.i) / sum66[idx6], w[idx6])
      term3 <- +ispx[, idx6, drop=FALSE] %*% Bv
      mat6 <- term1 + term2 + term3
    }

    ddpre[1:sk.i,1] <- matrix(t(mat1 + mat2 + mat3 + mat4 + mat5 + mat6), sk.i, 1)

    vec4 <- matrix(0, 1, k.i)
    if (length(idx4)) {
      A <- t(ispu[, idx4, drop=FALSE]) / sum44[idx4]
      A <- wrow(A, w[idx4])
      vec4 <- matrix(1,1,length(idx4)) %*% A
    }

    vec5 <- matrix(0, 1, k.i)
    if (length(idx5)) {
      A <- (t(ispv[, idx5, drop=FALSE]) - t(ispu[, idx5, drop=FALSE])) / sum55[idx5]
      A <- wrow(A, w[idx5])
      vec5 <- matrix(1,1,length(idx5)) %*% A
    }

    vec6 <- matrix(0, 1, k.i)
    if (length(idx6)) {
      A <- t(ispv[, idx6, drop=FALSE]) / sum66[idx6]
      A <- wrow(A, w[idx6])
      vec6 <- -matrix(1,1,length(idx6)) %*% A
    }

    ddpre[(sk.i+1):(sk.i+k.i), 1] <- vec4 + vec5 + vec6

    vec3 <- matrix(0, k.i, 1)
    if (length(idx3)) {
      A <- matrix(1, length(idx3), 1) / sum33[idx3]
      A <- wrow(A, w[idx3])
      vec3 <- mspx[, idx3, drop=FALSE] %*% A
    }

    vec6b <- matrix(0, k.i, 1)
    if (length(idx6)) {
      A <- matrix(1, length(idx6), 1) / sum66[idx6]
      A <- wrow(A, w[idx6])
      vec6b <- -ispx[, idx6, drop=FALSE] %*% A
    }

    ddpre[(sk.i+k.i+1):(sk.i+2*k.i), 1] <- vec3 + vec6b
    ddpre
  }


  hh_pre_w <- function(data1,data2,data3,data4,data5,data6,
                       sum11,sum22,sum33,sum44,sum55,sum66,
                       mat11,mat22,mat33,mat44,mat55,mat66,
                       vec44,vec55,vec66,vec33,vec66new,
                       w)
  {
    hhpre <- matrix(0, (sk.i+2*k.i), (sk.i+2*k.i))

    idx1 <- which(data1 != 0); idx2 <- which(data2 != 0); idx3 <- which(data3 != 0)
    idx4 <- which(data4 != 0); idx5 <- which(data5 != 0); idx6 <- which(data6 != 0)

    wrow <- function(A, ww) sweep(A, 1, ww, `*`)

    mat111 <- matrix(0, sk.i, sk.i)
    if (length(idx1)) {
      A <- t(mat11) / (sum11[idx1]^2)
      A <- wrow(A, w[idx1])
      mat111 <- mat11 %*% A
    }

    mat222 <- matrix(0, sk.i, sk.i)
    if (length(idx2)) {
      A <- t(mat22) / (sum22[idx2]^2)
      A <- wrow(A, w[idx2])
      mat222 <- mat22 %*% A
    }

    mat333 <- matrix(0, sk.i, sk.i)
    if (length(idx3)) {
      A <- t(mat33) / (sum33[idx3]^2)
      A <- wrow(A, w[idx3])
      mat333 <- mat33 %*% A
    }

    mat444 <- matrix(0, sk.i, sk.i)
    if (length(idx4)) {
      A <- t(mat44) / (sum44[idx4]^2)
      A <- wrow(A, w[idx4])
      mat444 <- mat44 %*% A
    }

    mat555 <- matrix(0, sk.i, sk.i)
    if (length(idx5)) {
      A <- t(mat55) / (sum55[idx5]^2)
      A <- wrow(A, w[idx5])
      mat555 <- mat55 %*% A
    }

    mat666 <- matrix(0, sk.i, sk.i)
    if (length(idx6)) {
      A <- t(mat66) / (sum66[idx6]^2)
      A <- wrow(A, w[idx6])
      mat666 <- mat66 %*% A
    }

    hhpre[1:sk.i, 1:sk.i] <- -mat111 - mat222 - mat333 - mat444 - mat555 - mat666

    vec444 <- matrix(0, sk.i, k.i)
    if (length(idx4)) {
      A <- t(vec44) / (sum44[idx4]^2)
      A <- wrow(A, w[idx4])
      vec444 <- mat44 %*% A
    }

    vec555 <- matrix(0, sk.i, k.i)
    if (length(idx5)) {
      A <- t(vec55) / (sum55[idx5]^2)
      A <- wrow(A, w[idx5])
      vec555 <- mat55 %*% A
    }

    vec666 <- matrix(0, sk.i, k.i)
    if (length(idx6)) {
      A <- t(vec66) / (sum66[idx6]^2)
      A <- wrow(A, w[idx6])
      vec666 <- mat66 %*% A
    }

    hhpre[1:sk.i, (sk.i+1):(sk.i+k.i)] <- -vec444 - vec555 - vec666

    vec333 <- matrix(0, sk.i, k.i)
    if (length(idx3)) {
      A <- t(vec33) / (sum33[idx3]^2)
      A <- wrow(A, w[idx3])
      vec333 <- mat33 %*% A
    }

    vec666b <- matrix(0, sk.i, k.i)
    if (length(idx6)) {
      A <- t(vec66new) / (sum66[idx6]^2)
      A <- wrow(A, w[idx6])
      vec666b <- mat66 %*% A
    }

    hhpre[1:sk.i, (sk.i+k.i+1):(sk.i+2*k.i)] <- -vec333 - vec666b

    vec4442 <- matrix(0, k.i, k.i)
    if (length(idx4)) {
      A <- t(vec44) / (sum44[idx4]^2)
      A <- wrow(A, w[idx4])
      vec4442 <- vec44 %*% A
    }

    vec5552 <- matrix(0, k.i, k.i)
    if (length(idx5)) {
      A <- t(vec55) / (sum55[idx5]^2)
      A <- wrow(A, w[idx5])
      vec5552 <- vec55 %*% A
    }

    vec6662 <- matrix(0, k.i, k.i)
    if (length(idx6)) {
      A <- t(vec66) / (sum66[idx6]^2)
      A <- wrow(A, w[idx6])
      vec6662 <- vec66 %*% A
    }

    hhpre[(sk.i+1):(sk.i+k.i), (sk.i+1):(sk.i+k.i)] <- -vec4442 - vec5552 - vec6662

    vec6663 <- matrix(0, k.i, k.i)
    if (length(idx6)) {
      A <- t(vec66new) / (sum66[idx6]^2)
      A <- wrow(A, w[idx6])
      vec6663 <- vec66 %*% A
    }
    hhpre[(sk.i+1):(sk.i+k.i), (sk.i+k.i+1):(sk.i+2*k.i)] <- -vec6663

    vec3333 <- matrix(0, k.i, k.i)
    if (length(idx3)) {
      A <- t(vec33) / (sum33[idx3]^2)
      A <- wrow(A, w[idx3])
      vec3333 <- vec33 %*% A
    }

    vec6664 <- matrix(0, k.i, k.i)
    if (length(idx6)) {
      A <- t(vec66new) / (sum66[idx6]^2)
      A <- wrow(A, w[idx6])
      vec6664 <- vec66new %*% A
    }

    hhpre[(sk.i+k.i+1):(sk.i+2*k.i), (sk.i+k.i+1):(sk.i+2*k.i)] <- -vec3333 - vec6664

    hhpre[(sk.i+1):(sk.i+k.i), 1:sk.i] <- t(hhpre[1:sk.i, (sk.i+1):(sk.i+k.i)])
    hhpre[(sk.i+k.i+1):(sk.i+2*k.i), 1:sk.i] <- t(hhpre[1:sk.i, (sk.i+k.i+1):(sk.i+2*k.i)])
    hhpre[(sk.i+k.i+1):(sk.i+2*k.i), (sk.i+1):(sk.i+k.i)] <- t(hhpre[(sk.i+1):(sk.i+k.i), (sk.i+k.i+1):(sk.i+2*k.i)])
    hhpre
  }


  get_knots <- function(x1, ctu, ctv)
  {
    probs <- seq_len(k.i - 3L) / (k.i - 2L)

    q1 <- as.numeric(
      quantile(x1, probs = probs, names = FALSE)
    )

    ## ctu=ctv=0 is used by the ORIGINAL terminal-event recoding
    ## when follow-up ends before the first secondary assessment.
    ## Zero is the boundary, not an interior inspection time, so exclude
    ## these artificial zeros when choosing secondary interior knots.
    knot_pre_2 <- c(ctu, ctv)
    knot_pre_2 <- knot_pre_2[is.finite(knot_pre_2) & knot_pre_2 > 0]

    if (length(knot_pre_2) < (k.i - 3L)) {
      stop("Too few positive secondary inspection endpoints for knot placement")
    }

    q2 <- as.numeric(
      quantile(knot_pre_2, probs = probs, names = FALSE)
    )

    ## Cover every finite likelihood-evaluation time and tau.
    qmax1 <- max(c(x1, tau), na.rm = TRUE) + 0.5
    qmax2 <- max(c(ctu, ctv, tau), na.rm = TRUE) + 0.5

    knot1 <- c(
      rep(0, spline.ord),
      q1,
      rep(qmax1, spline.ord)
    )

    knot2 <- c(
      rep(0, spline.ord),
      q2,
      rep(qmax2, spline.ord)
    )

    stopifnot(
      length(knot1) == n.spline + spline.ord + 1L,
      length(knot2) == n.spline + spline.ord + 1L
    )

    list(knot1 = knot1, knot2 = knot2)
  }


  totest_fast_diag <- function(x1, ctu, ctv,
                          delta1_1, delta2_1,
                          delta1_2, delta2_2, delta3_2,
                          w = rep(1, length(x1)),
                          alpha_init = NULL,
                          max_halvings = 10L)
  {
    knot12 <- get_knots(x1, ctu, ctv)
    knot1 <- knot12$knot1
    knot2 <- knot12$knot2

    ## I-spline / M-spline for primary margin.
    knotb <- knot1
    ord <- spline.ord

    bspx <- t(
      splines::splineDesign(
        knots = knotb,
        x = x1,
        ord = ord,
        outer.ok = TRUE
      )
    )

    pre.ispx <- apply(bspx, 2, rev)
    ispx <- apply(pre.ispx, 2, cumsum)
    ispx <- apply(ispx, 2, rev)
    ispx <- ispx[-1, , drop = FALSE]

    knotbM <- knot1[-c(1, n.spline + spline.ord + 1L)]
    ordM <- spline.ord - 1L

    bspxM <- t(
      splines::splineDesign(
        knots = knotbM,
        x = x1,
        ord = ordM,
        outer.ok = TRUE
      )
    )

    mspx <- bspxM /
      (knotbM[(1 + ordM):(n.spline + ordM)] -
         knotbM[1:n.spline]) * ordM

    ## I-spline for secondary margin at lower/upper endpoints.
    knotb <- knot2
    ord <- spline.ord

    bspu <- t(
      splines::splineDesign(
        knots = knotb,
        x = ctu,
        ord = ord,
        outer.ok = TRUE
      )
    )

    pre.ispu <- apply(bspu, 2, rev)
    ispu <- apply(pre.ispu, 2, cumsum)
    ispu <- apply(ispu, 2, rev)
    ispu <- ispu[-1, , drop = FALSE]

    bspv <- t(
      splines::splineDesign(
        knots = knotb,
        x = ctv,
        ord = ord,
        outer.ok = TRUE
      )
    )

    pre.ispv <- apply(bspv, 2, rev)
    ispv <- apply(pre.ispv, 2, cumsum)
    ispv <- apply(ispv, 2, rev)
    ispv <- ispv[-1, , drop = FALSE]

    ## ORIGINAL six observed-data categories.
    data1 <- delta1_1 * delta1_2
    data2 <- delta1_1 * delta2_2
    data3 <- delta1_1 * delta3_2
    data4 <- delta2_1 * delta1_2
    data5 <- delta2_1 * delta2_2
    data6 <- delta2_1 * delta3_2

    noze1 <- sum(data1)
    noze2 <- sum(data2)
    noze3 <- sum(data3)
    noze4 <- sum(data4)
    noze5 <- sum(data5)
    noze6 <- sum(data6)

    ## ORIGINAL derivative-design arrays.
    mat11 <- matrix(0, 1, noze1)
    mat22 <- matrix(0, 1, noze2)
    mat33 <- matrix(0, 1, noze3)
    mat44 <- matrix(0, 1, noze4)
    mat55 <- matrix(0, 1, noze5)
    mat66 <- matrix(0, 1, noze6)

    for (m in seq_len(k.i))
    {
      mat11 <- rbind(
        mat11,
        t(c(mspx[m, data1 != 0]) *
            t(ispu[, data1 != 0, drop = FALSE]))
      )

      mat22 <- rbind(
        mat22,
        t(c(mspx[m, data2 != 0]) *
            t(ispv[, data2 != 0, drop = FALSE])) -
          t(c(mspx[m, data2 != 0]) *
              t(ispu[, data2 != 0, drop = FALSE]))
      )

      mat33 <- rbind(
        mat33,
        t(c(mspx[m, data3 != 0]) *
            matrix(1, noze3, k.i)) -
          t(c(mspx[m, data3 != 0]) *
              t(ispv[, data3 != 0, drop = FALSE]))
      )

      mat44 <- rbind(
        mat44,
        ispu[, data4 != 0, drop = FALSE] -
          t(c(ispx[m, data4 != 0]) *
              t(ispu[, data4 != 0, drop = FALSE]))
      )

      mat55 <- rbind(
        mat55,
        ispv[, data5 != 0, drop = FALSE] -
          ispu[, data5 != 0, drop = FALSE] -
          t(c(ispx[m, data5 != 0]) *
              t(ispv[, data5 != 0, drop = FALSE])) +
          t(c(ispx[m, data5 != 0]) *
              t(ispu[, data5 != 0, drop = FALSE]))
      )

      mat66 <- rbind(
        mat66,
        -ispv[, data6 != 0, drop = FALSE] -
          t(c(ispx[m, data6 != 0]) *
              matrix(1, noze6, k.i)) +
          t(c(ispx[m, data6 != 0]) *
              t(ispv[, data6 != 0, drop = FALSE]))
      )
    }

    mat11 <- mat11[-1, , drop = FALSE]
    mat22 <- mat22[-1, , drop = FALSE]
    mat33 <- mat33[-1, , drop = FALSE]
    mat44 <- mat44[-1, , drop = FALSE]
    mat55 <- mat55[-1, , drop = FALSE]
    mat66 <- mat66[-1, , drop = FALSE]

    vec44 <- ispu[, data4 != 0, drop = FALSE]
    vec55 <- ispv[, data5 != 0, drop = FALSE] -
      ispu[, data5 != 0, drop = FALSE]
    vec66 <- -ispv[, data6 != 0, drop = FALSE]
    vec33 <- mspx[, data3 != 0, drop = FALSE]
    vec66new <- -ispx[, data6 != 0, drop = FALSE]

    ## Start.
    if (is.null(alpha_init)) {
      alpha <- matrix(0.5 / p_alpha, p_alpha, 1)
    } else {
      avec <- as.numeric(alpha_init)

      if (length(avec) != p_alpha || any(!is.finite(avec))) {
        stop("Invalid alpha_init")
      }

      alpha <- matrix(avec, p_alpha, 1)
      alpha[alpha < 0] <- 0

      if (sum(alpha) >= 1) {
        alpha <- alpha / sum(alpha) * (1 - 1e-10)
      }
    }

    current_parts <- sumvectorfunction(
      alpha, ispx, mspx, ispu, ispv
    )

    sum11 <- current_parts$sum11
    sum22 <- current_parts$sum22
    sum33 <- current_parts$sum33
    sum44 <- current_parts$sum44
    sum55 <- current_parts$sum55
    sum66 <- current_parts$sum66

    summl1 <- sumfunction_w(
      sum11, sum22, sum33, sum44, sum55, sum66,
      data1, data2, data3, data4, data5, data6,
      w
    )

    if (!is.finite(summl1)) {
      return(list(
        alpha = -9999,
        knot1 = -9999,
        knot2 = -9999,
        status = "invalid_start_likelihood",
        detail = "Initial observed-cell probability nonpositive/nonfinite"
      ))
    }

    dd <- gg_pre_w(
      ispx, mspx, ispu, ispv,
      data1, data2, data3, data4, data5, data6,
      noze1, noze2, noze3, noze4, noze5, noze6,
      sum11, sum22, sum33, sum44, sum55, sum66,
      w
    )

    hhh <- hh_pre_w(
      data1, data2, data3, data4, data5, data6,
      sum11, sum22, sum33, sum44, sum55, sum66,
      mat11, mat22, mat33, mat44, mat55, mat66,
      vec44, vec55, vec66, vec33, vec66new,
      w
    )

    ## Generalized gradient-projection algorithm:
    ## same statistical optimizer as the original code, with
    ## conditioning / active-set / line-search safeguards.
    indexA <- integer(0)
    aa <- NULL

    row_old <- -1L
    row_now <- 0L
    first_pass <- TRUE

    stopflag <- FALSE
    fail_reason <- "success"
    fail_detail <- NA_character_

    while (row_old != row_now)
    {
      if (first_pass)
      {
        first_pass <- FALSE

        hreg <- -hhh + diag(0.01, p_alpha)
        rc <- rcond(hreg)

        if (!is.finite(rc) || rc < 1e-15) {
          fail_reason <- "hessian_ill_conditioned"
          fail_detail <- paste0("initial rcond=", rc)
          stopflag <- TRUE
          break
        }

        ww <- solve(hreg)
        gt <- ww %*% dd

        row_now <- 0L
        pw <- diag(1, p_alpha)
        ddd <- pw %*% gt
        ddd[abs(ddd) < 2e-8] <- 0

      } else {

        if (row_now > 0L) {
          amat <- aa %*% ww %*% t(aa)
          rc <- rcond(amat)

          if (!is.finite(rc) || rc < 1e-15) {
            fail_reason <- "active_set_ill_conditioned"
            fail_detail <- paste0("rcond=", rc)
            stopflag <- TRUE
            break
          }

          pw <- diag(1, p_alpha) -
            ww %*% t(aa) %*% solve(amat) %*% aa

        } else {
          pw <- diag(1, p_alpha)
        }

        ddd <- pw %*% gt
        ddd[abs(ddd) < 2e-8] <- 0
      }

      inner_iter <- 0L

      while (max(abs(ddd)) > 1e-4)
      {
        inner_iter <- inner_iter + 1L

        if (inner_iter > 10000L) {
          fail_reason <- "iteration_limit"
          fail_detail <- "Exceeded 10000 inner iterations"
          stopflag <- TRUE
          break
        }

        tep <- numeric(0)
        tep1 <- integer(0)

        sumc <- sum(alpha)
        sumd <- sum(ddd)

        if ((1 - sumc) > 1e-7 && sumd > 1e-7) {
          tep <- c(tep, (1 - sumc) / sumd)
          tep1 <- c(tep1, -1L)
        }

        ## Joint coefficients use original 1e-7 numerical boundary.
        for (ii in seq_len(sk.i)) {
          if (alpha[ii] > 1e-7 && ddd[ii] < 0) {
            tep <- c(tep, (1e-7 - alpha[ii]) / ddd[ii])
            tep1 <- c(tep1, ii)
          }
        }

        ## Marginal coefficients use zero boundary.
        for (ii in (sk.i + 1L):p_alpha) {
          if (alpha[ii] > 0 && ddd[ii] < 0) {
            tep <- c(tep, -alpha[ii] / ddd[ii])
            tep1 <- c(tep1, ii)
          }
        }

        if (length(tep) > 0L) {
          kkkk <- which.min(tep)
          raw_step <- tep[kkkk]
          boundary_id <- tep1[kkkk]
        } else {
          kkkk <- NA_integer_
          raw_step <- Inf
          boundary_id <- NA_integer_
        }

        rrtt <- min(raw_step, 0.5)
        hit_boundary <- is.finite(raw_step) && raw_step <= 0.5

        summl2 <- summl1 - 1
        kk <- -1L
        alphat <- alpha

        while (summl2 <= summl1 && kk <= max_halvings)
        {
          kk <- kk + 1L
          alphat <- alpha + (0.5^kk) * rrtt * ddd

          if (sum(alphat) > 1 + 1e-5) {
            ## Invalid trial step: reject it and continue halving.
            ## This is a line-search event, not a dataset-level fit failure.
            summl2 <- -Inf
            next
          }

          add <- integer(0)

          if (hit_boundary && kk == 0L)
          {
            if (boundary_id != -1L) {
              if (boundary_id > sk.i) {
                alphat[boundary_id] <- 0
              } else {
                alphat[boundary_id] <- 1e-7
              }
            }

            for (ii in seq_len(p_alpha))
            {
              near_boundary <-
                (ii <= sk.i &&
                   abs(alphat[ii] - 1e-7) < 1e-8) ||
                (ii > sk.i && alphat[ii] < 1e-8)

              if (near_boundary &&
                  ii != boundary_id &&
                  !(ii %in% indexA)) {
                add <- c(add, ii)
              }
            }
          }

          alphat[alphat < 1e-8] <- 0

          trial_parts <- sumvectorfunction(
            alphat, ispx, mspx, ispu, ispv
          )

          summl2 <- sumfunction_w(
            trial_parts$sum11,
            trial_parts$sum22,
            trial_parts$sum33,
            trial_parts$sum44,
            trial_parts$sum55,
            trial_parts$sum66,
            data1, data2, data3, data4, data5, data6,
            w
          )

          if (is.infinite(summl2) && summl2 > 0) {
            fail_reason <- "trial_loglik_posInf"
            fail_detail <- "Trial log-likelihood was +Inf"
            stopflag <- TRUE
            break
          }
        }

        if (stopflag) break

        if (kk > max_halvings) {
          fail_reason <- "line_search_failure"
          fail_detail <- paste0(
            "No likelihood increase after >",
            max_halvings,
            " step halvings"
          )
          stopflag <- TRUE
          break
        }

        if (hit_boundary && kk == 0L)
        {
          if (boundary_id == -1L) {
            if (!(-1L %in% indexA)) {
              anew <- matrix(1, 1, p_alpha)
              indexA <- c(indexA, -1L)
              aa <- if (is.null(aa)) anew else rbind(aa, anew)
            }
          } else {
            if (!(boundary_id %in% indexA)) {
              anew <- matrix(0, 1, p_alpha)
              anew[1, boundary_id] <- -1
              indexA <- c(indexA, boundary_id)
              aa <- if (is.null(aa)) anew else rbind(aa, anew)
            }
          }

          if (length(add) > 0L) {
            for (jj in add) {
              if (!(jj %in% indexA)) {
                anew <- matrix(0, 1, p_alpha)
                anew[1, jj] <- -1
                indexA <- c(indexA, jj)
                aa <- if (is.null(aa)) anew else rbind(aa, anew)
              }
            }
          }

          row_now <- length(indexA)
        }

        alpha <- alphat

        current_parts <- sumvectorfunction(
          alpha, ispx, mspx, ispu, ispv
        )

        sum11 <- current_parts$sum11
        sum22 <- current_parts$sum22
        sum33 <- current_parts$sum33
        sum44 <- current_parts$sum44
        sum55 <- current_parts$sum55
        sum66 <- current_parts$sum66

        summl1 <- sumfunction_w(
          sum11, sum22, sum33, sum44, sum55, sum66,
          data1, data2, data3, data4, data5, data6,
          w
        )

        if (!is.finite(summl1)) {
          fail_reason <- "invalid_updated_likelihood"
          fail_detail <- "Accepted point has invalid likelihood"
          stopflag <- TRUE
          break
        }

        dd <- gg_pre_w(
          ispx, mspx, ispu, ispv,
          data1, data2, data3, data4, data5, data6,
          noze1, noze2, noze3, noze4, noze5, noze6,
          sum11, sum22, sum33, sum44, sum55, sum66,
          w
        )

        hhh <- hh_pre_w(
          data1, data2, data3, data4, data5, data6,
          sum11, sum22, sum33, sum44, sum55, sum66,
          mat11, mat22, mat33, mat44, mat55, mat66,
          vec44, vec55, vec66, vec33, vec66new,
          w
        )

        hreg <- -hhh + diag(0.01, p_alpha)
        rc <- rcond(hreg)

        if (!is.finite(rc) || rc < 1e-15) {
          fail_reason <- "hessian_ill_conditioned"
          fail_detail <- paste0("rcond=", rc)
          stopflag <- TRUE
          break
        }

        ww <- solve(hreg)
        gt <- ww %*% dd

        if (row_now > 0L) {
          amat <- aa %*% ww %*% t(aa)
          rc <- rcond(amat)

          if (!is.finite(rc) || rc < 1e-15) {
            fail_reason <- "active_set_ill_conditioned"
            fail_detail <- paste0("rcond=", rc)
            stopflag <- TRUE
            break
          }

          pw <- diag(1, p_alpha) -
            ww %*% t(aa) %*% solve(amat) %*% aa
        } else {
          pw <- diag(1, p_alpha)
        }

        ddd <- pw %*% gt
        ddd[abs(ddd) < 2e-8] <- 0
      }

      if (stopflag) break

      row_old <- row_now

      if (row_now > 0L)
      {
        amat <- aa %*% ww %*% t(aa)
        rc <- rcond(amat)

        if (!is.finite(rc) || rc < 1e-15) {
          fail_reason <- "active_set_ill_conditioned"
          fail_detail <- paste0("rcond=", rc)
          stopflag <- TRUE
          break
        }

        lambda <- solve(amat) %*% aa %*% gt

        if (min(lambda) < 0)
        {
          ik <- which.min(lambda)

          indexA <- indexA[-ik]

          if (length(indexA) == 0L) {
            aa <- NULL
          } else {
            aa <- aa[-ik, , drop = FALSE]
          }

          row_now <- length(indexA)
        }
      }
    }

    if (!stopflag) {
      list(
        alpha = alpha,
        knot1 = knot1,
        knot2 = knot2,
        status = "success",
        detail = NA_character_
      )
    } else {
      list(
        alpha = -9999,
        knot1 = -9999,
        knot2 = -9999,
        status = fail_reason,
        detail = fail_detail
      )
    }
  }


  validate_fit <- function(ans, tol = 1e-7)
  {
    if (is.null(ans$status)) {
      ans$status <- "invalid_fit_object"
      ans$detail <- "Missing status"
      return(ans)
    }

    if (!identical(ans$status, "success")) return(ans)

    a <- as.numeric(ans$alpha)

    if (length(a) != p_alpha || any(!is.finite(a))) {
      ans$status <- "invalid_coefficient_vector"
      ans$detail <- "Wrong length/nonfinite coefficient"
      return(ans)
    }

    if (any(a < -tol)) {
      ans$status <- "nonnegativity_constraint_failure"
      ans$detail <- paste0("min(alpha)=", min(a))
      return(ans)
    }

    if (sum(a) > 1 + tol) {
      ans$status <- "mass_constraint_failure"
      ans$detail <- paste0("sum(alpha)=", sum(a))
      return(ans)
    }

    if (any(!is.finite(ans$knot1)) ||
        any(!is.finite(ans$knot2))) {
      ans$status <- "invalid_knot_vector"
      ans$detail <- "Nonfinite knot"
      return(ans)
    }

    ans
  }


  softmax_alpha_from_eta <- function(eta)
  {
    eta <- as.numeric(eta)

    if (length(eta) != p_alpha || any(!is.finite(eta))) {
      stop("Invalid eta")
    }

    lower <- c(
      rep(1e-7, sk.i),
      rep(0, 2L*k.i)
    )

    R <- 1 - sum(lower)

    ## Stable softmax over c(eta, 0), where the last category is unused mass.
    shift <- max(c(eta, 0))
    ee <- exp(eta - shift)
    e0 <- exp(-shift)
    den <- sum(ee) + e0

    q <- ee / den
    q0 <- e0 / den

    alpha <- lower + R * q

    list(
      alpha = alpha,
      q = q,
      q0 = q0,
      lower = lower,
      R = R
    )
  }


  softmax_eta_from_alpha <- function(alpha_init = NULL,
                                     interior_mix = 1e-6)
  {
    lower <- c(
      rep(1e-7, sk.i),
      rep(0, 2L*k.i)
    )

    R <- 1 - sum(lower)

    if (is.null(alpha_init)) {
      ## Generic strictly interior alpha start with total mass 0.50.
      a <- rep(0.5 / p_alpha, p_alpha)
    } else {
      a <- as.numeric(alpha_init)

      if (length(a) != p_alpha || any(!is.finite(a))) {
        stop("Invalid alpha_init")
      }

      a <- pmax(a, lower)
    }

    ## If numerical rounding puts the start too close to / beyond total mass 1,
    ## retain its direction but restore comfortable slack.
    if (sum(a) >= 0.999999) {
      free <- pmax(a - lower, 0)
      target_free <- 0.95 - sum(lower)

      if (sum(free) > 0) {
        a <- lower + free / sum(free) * target_free
      } else {
        a <- lower + target_free / p_alpha
      }
    }

    free <- pmax(a - lower, 0)
    slack <- 1 - sum(a)

    q_full <- c(free / R, slack / R)

    ## Protect the inverse logit from exact zeros.  This changes only the
    ## numerical starting value, not the likelihood being optimized.
    eps <- interior_mix
    q_full <- (1 - eps) * q_full +
      eps * rep(1 / (p_alpha + 1L), p_alpha + 1L)

    q_full <- q_full / sum(q_full)

    log(q_full[seq_len(p_alpha)] / q_full[p_alpha + 1L])
  }


  totest_bfgs_diag <- function(x1, ctu, ctv,
                               delta1_1, delta2_1,
                               delta1_2, delta2_2, delta3_2,
                               w = rep(1, length(x1)),
                               alpha_init = NULL,
                               maxit = 300L,
                               reltol = 1e-8)
  {
    knot12 <- get_knots(x1, ctu, ctv)
    knot1 <- knot12$knot1
    knot2 <- knot12$knot2

    ## Same I-spline / M-spline basis construction as the original estimator.
    knotb <- knot1
    ord <- spline.ord

    bspx <- t(
      splines::splineDesign(
        knots = knotb,
        x = x1,
        ord = ord,
        outer.ok = TRUE
      )
    )

    pre.ispx <- apply(bspx, 2, rev)
    ispx <- apply(pre.ispx, 2, cumsum)
    ispx <- apply(ispx, 2, rev)
    ispx <- ispx[-1, , drop = FALSE]

    knotbM <- knot1[-c(1, n.spline + spline.ord + 1L)]
    ordM <- spline.ord - 1L

    bspxM <- t(
      splines::splineDesign(
        knots = knotbM,
        x = x1,
        ord = ordM,
        outer.ok = TRUE
      )
    )

    mspx <- bspxM /
      (knotbM[(1 + ordM):(n.spline + ordM)] -
         knotbM[1:n.spline]) * ordM

    knotb <- knot2
    ord <- spline.ord

    bspu <- t(
      splines::splineDesign(
        knots = knotb,
        x = ctu,
        ord = ord,
        outer.ok = TRUE
      )
    )

    pre.ispu <- apply(bspu, 2, rev)
    ispu <- apply(pre.ispu, 2, cumsum)
    ispu <- apply(ispu, 2, rev)
    ispu <- ispu[-1, , drop = FALSE]

    bspv <- t(
      splines::splineDesign(
        knots = knotb,
        x = ctv,
        ord = ord,
        outer.ok = TRUE
      )
    )

    pre.ispv <- apply(bspv, 2, rev)
    ispv <- apply(pre.ispv, 2, cumsum)
    ispv <- apply(ispv, 2, rev)
    ispv <- ispv[-1, , drop = FALSE]

    ## Original six mixed-censoring cells.
    data1 <- delta1_1 * delta1_2
    data2 <- delta1_1 * delta2_2
    data3 <- delta1_1 * delta3_2
    data4 <- delta2_1 * delta1_2
    data5 <- delta2_1 * delta2_2
    data6 <- delta2_1 * delta3_2

    noze1 <- sum(data1)
    noze2 <- sum(data2)
    noze3 <- sum(data3)
    noze4 <- sum(data4)
    noze5 <- sum(data5)
    noze6 <- sum(data6)

    ## Cache because optim commonly asks for f and gr at the same eta.
    cache <- new.env(parent = emptyenv())
    cache$eta <- NULL
    cache$value <- NULL

    evaluate_eta <- function(eta)
    {
      eta <- as.numeric(eta)

      if (!is.null(cache$eta) &&
          length(cache$eta) == length(eta) &&
          isTRUE(all.equal(cache$eta, eta, tolerance = 0))) {
        return(cache$value)
      }

      tr <- softmax_alpha_from_eta(eta)
      alpha <- tr$alpha

      sv <- sumvectorfunction(
        alpha,
        ispx, mspx, ispu, ispv
      )

      ll <- sumfunction_w(
        sv$sum11, sv$sum22, sv$sum33,
        sv$sum44, sv$sum55, sv$sum66,
        data1, data2, data3,
        data4, data5, data6,
        w
      )

      if (!is.finite(ll)) {
        ans <- list(
          ll = -Inf,
          grad_eta_loglik = rep(NA_real_, p_alpha),
          alpha = alpha
        )
        cache$eta <- eta
        cache$value <- ans
        return(ans)
      }

      score_alpha <- as.numeric(
        gg_pre_w(
          ispx, mspx, ispu, ispv,
          data1, data2, data3,
          data4, data5, data6,
          noze1, noze2, noze3,
          noze4, noze5, noze6,
          sv$sum11, sv$sum22, sv$sum33,
          sv$sum44, sv$sum55, sv$sum66,
          w
        )
      )

      ## Chain rule:
      ## d alpha_j / d eta_k = R q_j (I[j=k] - q_k)
      ## therefore
      ## d ll / d eta_k = R q_k {score_k - sum_j q_j score_j}.
      q <- tr$q
      weighted_score_mean <- sum(q * score_alpha)

      grad_eta_loglik <-
        tr$R * q * (score_alpha - weighted_score_mean)

      ans <- list(
        ll = ll,
        grad_eta_loglik = grad_eta_loglik,
        alpha = alpha
      )

      cache$eta <- eta
      cache$value <- ans
      ans
    }

    eta0 <- tryCatch(
      softmax_eta_from_alpha(
        alpha_init = alpha_init,
        interior_mix = 1e-6
      ),
      error = function(e) e
    )

    if (inherits(eta0, "error")) {
      return(list(
        alpha = -9999,
        knot1 = -9999,
        knot2 = -9999,
        status = "bfgs_invalid_start",
        detail = conditionMessage(eta0)
      ))
    }

    z0 <- evaluate_eta(eta0)

    if (!is.finite(z0$ll) ||
        any(!is.finite(z0$grad_eta_loglik))) {
      return(list(
        alpha = -9999,
        knot1 = -9999,
        knot2 = -9999,
        status = "bfgs_invalid_start_likelihood",
        detail = "Softmax/BFGS start has invalid likelihood or gradient"
      ))
    }

    fn <- function(eta)
    {
      z <- evaluate_eta(eta)

      if (!is.finite(z$ll)) {
        return(1e100)
      }

      -z$ll
    }

    gr <- function(eta)
    {
      z <- evaluate_eta(eta)

      if (!is.finite(z$ll) ||
          any(!is.finite(z$grad_eta_loglik))) {
        return(rep(0, p_alpha))
      }

      -z$grad_eta_loglik
    }

    fit <- tryCatch(
      optim(
        par = eta0,
        fn = fn,
        gr = gr,
        method = "BFGS",
        control = list(
          maxit = maxit,
          reltol = reltol
        )
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(list(
        alpha = -9999,
        knot1 = -9999,
        knot2 = -9999,
        status = "bfgs_error",
        detail = conditionMessage(fit)
      ))
    }

    final <- evaluate_eta(fit$par)
    alpha_hat <- as.numeric(final$alpha)

    if (any(!is.finite(alpha_hat)) ||
        !is.finite(final$ll) ||
        any(!is.finite(final$grad_eta_loglik))) {
      return(list(
        alpha = -9999,
        knot1 = -9999,
        knot2 = -9999,
        status = "bfgs_invalid_final_solution",
        detail = "Nonfinite final coefficient, likelihood, or transformed gradient"
      ))
    }

    ## Defensive feasibility check of the transformed solution.
    if (any(alpha_hat[seq_len(sk.i)] < 1e-7 - 1e-12) ||
        any(alpha_hat[(sk.i + 1L):p_alpha] < -1e-12) ||
        sum(alpha_hat) >= 1 + 1e-12) {
      return(list(
        alpha = -9999,
        knot1 = -9999,
        knot2 = -9999,
        status = "bfgs_transform_feasibility_failure",
        detail = paste0(
          "min joint=", min(alpha_hat[seq_len(sk.i)]),
          "; min marginal=",
          min(alpha_hat[(sk.i + 1L):p_alpha]),
          "; sum=", sum(alpha_hat)
        )
      ))
    }

    if (!is.null(fit$convergence) && fit$convergence != 0) {
      return(list(
        alpha = -9999,
        knot1 = -9999,
        knot2 = -9999,
        status = paste0("bfgs_convergence_", fit$convergence),
        detail = paste0(
          fit$message,
          "; max|grad_eta|=",
          signif(max(abs(final$grad_eta_loglik)), 5),
          "; fn/gr evals=",
          paste(fit$counts, collapse = "/")
        )
      ))
    }

    list(
      alpha = matrix(alpha_hat, p_alpha, 1),
      knot1 = knot1,
      knot2 = knot2,
      status = "success",
      detail = paste0(
        "softmax/BFGS; max|grad_eta|=",
        signif(max(abs(final$grad_eta_loglik)), 5),
        "; fn/gr evals=",
        paste(fit$counts, collapse = "/")
      )
    )
  }


  safe_fit_bfgs <- function(dat,
                            w = rep(1, size),
                            alpha_init = NULL,
                            maxit = 300L,
                            reltol = 1e-8)
  {
    warning_messages <- character(0)

    ans <- tryCatch(
      withCallingHandlers(
        totest_bfgs_diag(
          dat$x1, dat$ctu, dat$ctv,
          dat$delta1_1, dat$delta2_1,
          dat$delta1_2, dat$delta2_2, dat$delta3_2,
          w = w,
          alpha_init = alpha_init,
          maxit = maxit,
          reltol = reltol
        ),
        warning = function(wrn) {
          warning_messages <<- c(
            warning_messages,
            conditionMessage(wrn)
          )
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        list(
          alpha = -9999,
          knot1 = -9999,
          knot2 = -9999,
          status = "R_error",
          detail = conditionMessage(e)
        )
      }
    )

    ans <- validate_fit(ans)

    ans$n_warnings <- length(warning_messages)

    ans$warning_detail <-
      if (length(warning_messages) == 0L) {
        NA_character_
      } else {
        paste(unique(warning_messages), collapse = " | ")
      }

    ans
  }


  safe_fit_fast <- function(dat,
                       w = rep(1, size),
                       alpha_init = NULL,
                       max_halvings = 10L)
  {
    warning_messages <- character(0)

    ans <- tryCatch(
      withCallingHandlers(
        totest_fast_diag(
          dat$x1, dat$ctu, dat$ctv,
          dat$delta1_1, dat$delta2_1,
          dat$delta1_2, dat$delta2_2, dat$delta3_2,
          w = w,
          alpha_init = alpha_init,
          max_halvings = max_halvings
        ),
        warning = function(wrn) {
          warning_messages <<- c(
            warning_messages,
            conditionMessage(wrn)
          )
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        list(
          alpha = -9999,
          knot1 = -9999,
          knot2 = -9999,
          status = "R_error",
          detail = conditionMessage(e)
        )
      }
    )

    ans <- validate_fit(ans)

    ans$n_warnings <- length(warning_messages)

    ans$warning_detail <-
      if (length(warning_messages) == 0L) {
        NA_character_
      } else {
        paste(unique(warning_messages), collapse = " | ")
      }

    ans
  }


  make_interior_start <- function(alpha_hat, mix = 1e-4)
  {
    generic <- matrix(0.5 / p_alpha, p_alpha, 1)

    ah <- matrix(as.numeric(alpha_hat), p_alpha, 1)
    ah[ah < 0] <- 0

    if (sum(ah) >= 1) {
      ah <- ah / sum(ah) * (1 - 1e-10)
    }

    (1 - mix) * ah + mix * generic
  }


  make_generic_start <- function(total_mass = 0.25)
  {
    stopifnot(total_mass > 0, total_mass < 1)
    matrix(total_mass / p_alpha, p_alpha, 1)
  }


  fit_original_hybrid <- function(dat)
  {
    ## Stage 1: one fast original active-set fit.
    fit1 <- safe_fit_fast(
      dat,
      w = rep(1, size),
      alpha_init = NULL,
      max_halvings = 30L
    )

    if (identical(fit1$status, "success")) {
      return(list(
        fit = fit1,
        initial_status = "success",
        initial_detail = NA_character_,
        fast_rescue_used = FALSE,
        fast_rescue_status = NA_character_,
        fast_rescue_detail = NA_character_,
        bfgs_fallback_used = FALSE,
        bfgs_fallback_status = NA_character_,
        bfgs_fallback_detail = NA_character_
      ))
    }

    ## Stage 2: softmax/BFGS fallback on the SAME generated dataset.
    fit2 <- safe_fit_bfgs(
      dat,
      w = rep(1, size),
      alpha_init = make_generic_start(total_mass = 0.25),
      maxit = 500L,
      reltol = 1e-8
    )

    list(
      fit = fit2,
      initial_status = fit1$status,
      initial_detail = fit1$detail,
      fast_rescue_used = FALSE,
      fast_rescue_status = NA_character_,
      fast_rescue_detail = NA_character_,
      bfgs_fallback_used = TRUE,
      bfgs_fallback_status = fit2$status,
      bfgs_fallback_detail = fit2$detail
    )
  }


  fit_bootstrap_with_rescue <- function(dat, w, alpha_hat)
  {
    ## Stage 1: one fast warm-start fit.
    fit1 <- safe_fit_fast(
      dat,
      w = w,
      alpha_init = make_interior_start(alpha_hat, mix = 1e-4),
      max_halvings = 30L
    )

    if (identical(fit1$status, "success")) {
      return(list(
        fit = fit1,
        initial_status = "success",
        initial_detail = NA_character_,
        fast_rescue_used = FALSE,
        fast_rescue_status = NA_character_,
        fast_rescue_detail = NA_character_,
        after_fast_status = "success",
        bfgs_fallback_used = FALSE,
        bfgs_fallback_status = NA_character_,
        bfgs_fallback_detail = NA_character_
      ))
    }

    ## Stage 2: SAME multiplier realization, softmax/BFGS fallback.
    ## This is NOT a replacement draw.
    fit2 <- safe_fit_bfgs(
      dat,
      w = w,
      alpha_init = make_interior_start(alpha_hat, mix = 0.10),
      maxit = 300L,
      reltol = 1e-8
    )

    list(
      fit = fit2,
      initial_status = fit1$status,
      initial_detail = fit1$detail,
      fast_rescue_used = FALSE,
      fast_rescue_status = NA_character_,
      fast_rescue_detail = NA_character_,
      after_fast_status = fit1$status,
      bfgs_fallback_used = TRUE,
      bfgs_fallback_status = fit2$status,
      bfgs_fallback_detail = fit2$detail
    )
  }


  ispline_eval <- function(t, knot_vec)
  {
    bsp <- t(splines::splineDesign(knots = knot_vec, x = t, ord = spline.ord, outer.ok = TRUE))
    rev(cumsum(rev(bsp)))[-1]
  }


  mspline_eval <- function(t, knot_vec)
  {
    knotb <- knot_vec[-c(1, (n.spline + spline.ord + 1))]
    ord <- spline.ord - 1
    bsp <- t(splines::splineDesign(knots = knotb, x = t, ord = ord, outer.ok = TRUE))
    as.numeric(bsp / (knotb[(1 + ord):(n.spline + ord)] - knotb[1:n.spline]) * ord)
  }


  joint_fit <- function(t1, t2, fit) {
    isp1 <- ispline_eval(t1, fit$knot1)
    isp2 <- ispline_eval(t2, fit$knot2)
    mu   <- t(matrix(fit$alpha[1:sk.i], k.i, k.i))
    drop(matrix(isp1, 1, n.spline) %*% mu %*% matrix(isp2, n.spline, 1))
  }


  mar1_fit <- function(t, fit) {
    isp <- ispline_eval(t, fit$knot1)
    mu  <- t(matrix(fit$alpha[1:sk.i], k.i, k.i))
    omg <- matrix(fit$alpha[(sk.i+k.i+1):(sk.i+2*k.i)], k.i, 1)
    drop(matrix(isp, 1, n.spline) %*% mu %*% matrix(1, n.spline, 1) +
           matrix(isp, 1, n.spline) %*% omg)
  }


  mar2_fit <- function(s, fit) {
    isp <- ispline_eval(s, fit$knot2)
    mu  <- t(matrix(fit$alpha[1:sk.i], k.i, k.i))
    pi  <- matrix(fit$alpha[(sk.i+1):(sk.i+k.i)], k.i, 1)
    drop(matrix(1, 1, n.spline) %*% mu %*% matrix(isp, n.spline, 1) +
           matrix(isp, 1, n.spline) %*% pi)
  }


  barF_fit <- function(t, s, fit) {
    val <- 1 - mar1_fit(t, fit) - mar2_fit(s, fit) + joint_fit(t, s, fit)
    max(min(val, 1), 0)
  }


  dF2_ds_fit <- function(s, fit) {
    msp2 <- mspline_eval(s, fit$knot2)
    mu   <- t(matrix(fit$alpha[1:sk.i], k.i, k.i))
    pi   <- matrix(fit$alpha[(sk.i+1):(sk.i+k.i)], k.i, 1)
    drop(matrix(1, 1, n.spline) %*% mu %*% matrix(msp2, n.spline, 1) +
           matrix(msp2, 1, n.spline) %*% pi)
  }


  dF_ds_fit <- function(t, s, fit) {
    isp1 <- ispline_eval(t, fit$knot1)
    msp2 <- mspline_eval(s, fit$knot2)
    mu   <- t(matrix(fit$alpha[1:sk.i], k.i, k.i))
    drop(matrix(isp1, 1, n.spline) %*% mu %*% matrix(msp2, n.spline, 1))
  }


  dbarF_ds_fit <- function(t, s, fit) {
    -(dF2_ds_fit(s, fit)) + dF_ds_fit(t, s, fit)
  }


  wr_integral_aim21_local <- function(fitC, fitT, tau,
                                      tau_l = 0,
                                      dt = 0.02, ds = 0.02) {
    stopifnot(is.finite(tau_l), tau_l >= 0, tau_l < tau,
              is.finite(dt), dt > 0, is.finite(ds), ds > 0)

    ## The exact/right-censored primary margin is identifiable from time zero.
    ## Preserve its ordinary finite-horizon comparison over [0, tau].
    t_edges <- seq(0, tau, length.out = ceiling(tau / dt) + 1L)
    t_right <- t_edges[-1L]
    F0_edges <- vapply(t_edges, mar1_fit, numeric(1), fit = fitC)
    F1_edges <- vapply(t_edges, mar1_fit, numeric(1), fit = fitT)
    dF0 <- diff(F0_edges)
    dF1 <- diff(F1_edges)
    ST1 <- pmax(0, pmin(1, 1 - vapply(t_right, mar1_fit,
                                     numeric(1), fit = fitT)))
    ST0 <- pmax(0, pmin(1, 1 - vapply(t_right, mar1_fit,
                                     numeric(1), fit = fitC)))
    W_primary <- sum(ST1 * dF0)
    L_primary <- sum(ST0 * dF1)

    ## The interval-censored secondary margin is restricted to [tau_l, tau].
    ## Midpoint quadrature avoids evaluation outside the restricted interval.
    s_edges <- seq(tau_l, tau,
                   length.out = ceiling((tau - tau_l) / ds) + 1L)
    s_mid <- (s_edges[-length(s_edges)] + s_edges[-1L]) / 2
    widths <- diff(s_edges)
    barF1_tau_s <- vapply(s_mid, function(s) barF_fit(tau, s, fitT), numeric(1))
    dbarF0_ds <- vapply(s_mid, function(s) dbarF_ds_fit(tau, s, fitC), numeric(1))
    barF0_tau_s <- vapply(s_mid, function(s) barF_fit(tau, s, fitC), numeric(1))
    dbarF1_ds <- vapply(s_mid, function(s) dbarF_ds_fit(tau, s, fitT), numeric(1))

    W_sec <- -sum(barF1_tau_s * dbarF0_ds * widths)
    L_sec <- -sum(barF0_tau_s * dbarF1_ds * widths)

    W <- W_primary + W_sec
    L <- L_primary + L_sec
    if (!all(is.finite(c(W, L, W_primary, W_sec, L_primary, L_sec))) ||
        W <= 0 || L <= 0 ||
        min(W_primary, W_sec, L_primary, L_sec) < -1e-7) {
      return(c(W=NA_real_, L=NA_real_, WR=NA_real_))
    }
    c(W=W, L=L, WR=W/L,
      W_primary=W_primary, W_sec=W_sec,
      L_primary=L_primary, L_sec=L_sec)
  }


  safe_wr_tau <- function(fitC, fitT, tau_l)
  {
    out <- tryCatch(
      wr_integral_aim21_local(
        fitC, fitT,
        tau = tau,
        tau_l = tau_l,
        dt = dt,
        ds = ds
      ),
      error = function(e) {
        c(W = NA_real_, L = NA_real_, WR = NA_real_)
      }
    )

    wr <- suppressWarnings(as.numeric(out["WR"]))
    if (!is.finite(wr) || wr <= 0) {
      return(list(status = "functional_failure", WR = NA_real_, logWR = NA_real_))
    }
    list(status = "success", WR = wr, logWR = log(wr))
  }

  configure <- function(n, k, tau_value) {
    size <<- as.integer(n)
    tau <<- as.numeric(tau_value)
    k.i <<- as.integer(k)
    n.spline <<- k.i
    sk.i <<- k.i^2
    spline.ord <<- 4L
    p_alpha <<- sk.i + 2L*k.i
    dt <<- 0.02
    ds <<- 0.02
    invisible(NULL)
  }

  fit_original <- function(dat, k, tau_value) {
    configure(length(dat$x1), k, tau_value)
    fit_original_hybrid(dat)
  }

  fit_bootstrap <- function(dat, w, alpha_hat, k, tau_value) {
    configure(length(dat$x1), k, tau_value)
    fit_bootstrap_with_rescue(dat, w, alpha_hat)
  }

  evaluate <- function(fit0, fit1, k, tau_value, tau_l = 0) {
    configure(1L, k, tau_value)
    out <- tryCatch(
      wr_integral_aim21_local(fit0, fit1, tau = tau, tau_l = tau_l, dt = dt, ds = ds),
      error = function(e) c(W=NA_real_, L=NA_real_, WR=NA_real_)
    )
    W <- suppressWarnings(as.numeric(out["W"]))
    L <- suppressWarnings(as.numeric(out["L"]))
    WR <- suppressWarnings(as.numeric(out["WR"]))
    if (!is.finite(W) || !is.finite(L) || !is.finite(WR) || W <= 0 || L <= 0 || WR <= 0) {
      return(list(status="functional_failure", detail=NA_character_, logwr=NA_real_, win=W, lose=L, wr=WR))
    }
    list(status="success", detail=NA_character_, logwr=log(WR), win=W, lose=L, wr=WR)
  }

  list(
    fit_original = fit_original,
    fit_bootstrap = fit_bootstrap,
    evaluate = evaluate
  )
})


#' Estimate a prioritized win ratio under interval and/or right censoring
#'
#' @param data A data.frame with one row per subject.
#' @param treatment Name of a binary treatment column. Numeric treatment must be 0/1.
#' @param margin1 List describing the higher-priority outcome. For type='interval',
#'   specify left/right censoring-bound columns; optional visit_lower/visit_upper
#'   columns retain the finite inspection pair used for spline-knot placement.
#' @param margin2 List describing the lower-priority outcome, with the same
#'   interval-margin convention as margin1.
#' @param tau Positive finite analysis horizon.
#' @param tau_l Lower integration bound. Must satisfy 0 <= tau_l < tau.
#'   For interval/interval, both prioritized win/loss integrals use [tau_l, tau].
#'   For right/interval, the primary comparison remains on [0, tau] and only
#'   the interval-censored secondary comparison uses [tau_l, tau]. tau_l=0
#'   computes the untrimmed exploratory target; the positive-lower-bound
#'   asymptotic theory applies when tau_l is strictly positive and satisfies
#'   the assumptions of the methodology.
#' @param B Number of fixed multiplier bootstrap draws.
#' @param k Fixed sieve dimension (default 6).
#' @param parallel Reserved for a future parallel implementation. Keep FALSE in this version.
#' @param seed Optional bootstrap random seed.
#' @param conf_level Confidence level for the log-scale Wald interval.
#' @return An object of class winratio_ic.
#' @author Zhiguo Li <zhiguo.li@duke.edu>
winratio_ic <- function(data,
                        treatment,
                        margin1,
                        margin2,
                        tau,
                        tau_l = 0,
                        B = 200L,
                        k = 6L,
                        parallel = FALSE,
                        seed = NULL,
                        conf_level = 0.95) {

  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (!nrow(data)) stop("data has no rows.", call. = FALSE)
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("tau must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(tau_l) || length(tau_l) != 1L || !is.finite(tau_l) ||
      tau_l < 0 || tau_l >= tau) {
    stop("tau_l must be one finite number satisfying 0 <= tau_l < tau.", call. = FALSE)
  }
  B <- as.integer(B)
  k <- as.integer(k)
  if (is.na(B) || B < 2L) stop("B must be an integer >= 2.", call. = FALSE)
  if (is.na(k) || k < 4L) stop("k must be an integer >= 4.", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be between 0 and 1.", call. = FALSE)
  }
  if (isTRUE(parallel)) {
    stop("parallel=TRUE is not implemented in this first companion-software version; use parallel=FALSE.", call. = FALSE)
  }

  type1 <- .wr_validate_margin(margin1, "margin1")
  type2 <- .wr_validate_margin(margin2, "margin2")

  scheme <- paste(type1, type2, sep = "/")
  if (!scheme %in% c("interval/interval", "right/interval")) {
    if (scheme == "interval/right") {
      stop(
        "The supplied mixed-censoring source code supports a right-censored primary margin and an interval-censored secondary margin. ",
        "The reverse interval/right likelihood is not implemented from the supplied source and is therefore not silently approximated.",
        call. = FALSE
      )
    }
    stop("Unsupported censoring combination: ", scheme, call. = FALSE)
  }

  groups <- .wr_split_treatment(data, treatment)
  if (!length(groups$control) || !length(groups$treatment)) {
    stop("Both treatment groups must contain at least one subject.", call. = FALSE)
  }

  if (scheme == "interval/interval") {
    m1 <- .wr_prepare_interval(data, margin1, "margin1")
    m2 <- .wr_prepare_interval(data, margin2, "margin2")

    make_arm <- function(idx) {
      a <- .wr_subset_prepared(m1, idx)
      b <- .wr_subset_prepared(m2, idx)
      list(
        tu_1 = a$tu, tu_2 = b$tu,
        tv_1 = a$tv, tv_2 = b$tv,
        delta1_1 = a$delta1, delta2_1 = a$delta2, delta3_1 = a$delta3,
        delta1_2 = b$delta1, delta2_2 = b$delta2, delta3_2 = b$delta3
      )
    }
    dat0 <- make_arm(groups$control)
    dat1 <- make_arm(groups$treatment)
    engine <- .wr_ii_engine

  } else {
    m1 <- .wr_prepare_right(data, margin1, "margin1")
    m2 <- .wr_prepare_interval(data, margin2, "margin2")

    make_arm <- function(idx) {
      a <- .wr_subset_prepared(m1, idx)
      b <- .wr_subset_prepared(m2, idx)
      list(
        x1 = a$time,
        ctu = b$tu,
        ctv = b$tv,
        delta1_1 = a$delta_event,
        delta2_1 = a$delta_censor,
        delta1_2 = b$delta1,
        delta2_2 = b$delta2,
        delta3_2 = b$delta3
      )
    }
    dat0 <- make_arm(groups$control)
    dat1 <- make_arm(groups$treatment)
    engine <- .wr_ri_engine
  }

  original0 <- engine$fit_original(dat0, k, tau)
  original1 <- engine$fit_original(dat1, k, tau)
  fit0 <- original0$fit
  fit1 <- original1$fit

  if (!identical(fit0$status, "success")) {
    stop("Control-arm fit failed: ", fit0$status, if (!is.null(fit0$detail) && !is.na(fit0$detail)) paste0(" (", fit0$detail, ")") else "", call. = FALSE)
  }
  if (!identical(fit1$status, "success")) {
    stop("Treatment-arm fit failed: ", fit1$status, if (!is.null(fit1$detail) && !is.na(fit1$detail)) paste0(" (", fit1$detail, ")") else "", call. = FALSE)
  }

  point <- engine$evaluate(fit0, fit1, k, tau, tau_l)
  if (!identical(point$status, "success") || !is.finite(point$logwr)) {
    stop("Win-ratio functional failed after successful arm fits: ", point$status, call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  n0 <- length(groups$control)
  n1 <- length(groups$treatment)
  W0 <- matrix(stats::rexp(B * n0, rate = 1), nrow = B, ncol = n0)
  W1 <- matrix(stats::rexp(B * n1, rate = 1), nrow = B, ncol = n1)
  W0 <- W0 / rowMeans(W0)
  W1 <- W1 / rowMeans(W1)

  boot_logwr <- rep(NA_real_, B)
  boot_wr <- rep(NA_real_, B)
  boot_status <- rep(NA_character_, B)
  fit0_status <- rep(NA_character_, B)
  fit1_status <- rep(NA_character_, B)
  rescue0 <- rep(FALSE, B)
  rescue1 <- rep(FALSE, B)

  for (b in seq_len(B)) {
    r0 <- engine$fit_bootstrap(dat0, W0[b, ], fit0$alpha, k, tau)
    r1 <- engine$fit_bootstrap(dat1, W1[b, ], fit1$alpha, k, tau)

    fb0 <- r0$fit
    fb1 <- r1$fit
    fit0_status[b] <- fb0$status
    fit1_status[b] <- fb1$status

    rescue0[b] <- isTRUE(r0$rescue_used) || isTRUE(r0$bfgs_fallback_used)
    rescue1[b] <- isTRUE(r1$rescue_used) || isTRUE(r1$bfgs_fallback_used)

    if (!identical(fb0$status, "success") || !identical(fb1$status, "success")) {
      boot_status[b] <- "fit_failure"
      next
    }

    z <- engine$evaluate(fb0, fb1, k, tau, tau_l)
    if (!identical(z$status, "success") || !is.finite(z$logwr)) {
      boot_status[b] <- "functional_failure"
      next
    }

    boot_status[b] <- "success"
    boot_logwr[b] <- z$logwr
    boot_wr[b] <- z$wr
  }

  ok <- boot_status == "success" & is.finite(boot_logwr)
  n_ok <- sum(ok)
  se_log <- if (n_ok >= 2L) stats::sd(boot_logwr[ok]) else NA_real_

  alpha <- 1 - conf_level
  zcrit <- stats::qnorm(1 - alpha/2)
  ci_log <- if (is.finite(se_log)) {
    point$logwr + c(-1, 1) * zcrit * se_log
  } else {
    c(NA_real_, NA_real_)
  }
  ci_wr <- exp(ci_log)

  boot_diag <- data.frame(
    bootstrap_rep = seq_len(B),
    status = boot_status,
    control_fit_status = fit0_status,
    treatment_fit_status = fit1_status,
    control_rescue_used = rescue0,
    treatment_rescue_used = rescue1,
    win_ratio = boot_wr,
    log_win_ratio = boot_logwr,
    stringsAsFactors = FALSE
  )

  ans <- list(
    call = match.call(),
    scheme = scheme,
    margin1_priority = "higher",
    control_label = groups$control_label,
    treatment_label = groups$treatment_label,
    n = c(control = n0, treatment = n1),
    tau = tau,
    tau_l = tau_l,
    estimator = if (tau_l > 0) "trimmed" else "untrimmed",
    k = k,
    B = B,
    win_probability = point$win,
    loss_probability = point$lose,
    win_ratio = point$wr,
    log_win_ratio = point$logwr,
    se_log_win_ratio = se_log,
    conf_level = conf_level,
    ci_log_win_ratio = setNames(ci_log, c("lower", "upper")),
    ci_win_ratio = setNames(ci_wr, c("lower", "upper")),
    bootstrap = list(
      n_attempted = B,
      n_success = n_ok,
      success_fraction = n_ok / B,
      diagnostics = boot_diag
    ),
    fit = list(control = fit0, treatment = fit1)
  )

  class(ans) <- "winratio_ic"
  ans
}

print.winratio_ic <- function(x, ...) {
  cat("Win-ratio analysis\n")
  cat("Censoring scheme:  ", x$scheme, "\n", sep = "")
  cat("Control / treatment n: ", x$n[["control"]], " / ", x$n[["treatment"]], "\n", sep = "")
  cat("Analysis horizon tau: ", format(x$tau), "\n", sep = "")
  cat("Lower bound tau_l:    ", format(x$tau_l), "\n", sep = "")
  cat("Estimator:            ", x$estimator, "\n", sep = "")
  cat("Win probability:    ", formatC(x$win_probability, digits = 4, format = "f"), "\n", sep = "")
  cat("Loss probability:   ", formatC(x$loss_probability, digits = 4, format = "f"), "\n", sep = "")
  cat("Win ratio:          ", formatC(x$win_ratio, digits = 4, format = "f"), "\n", sep = "")
  cat("log(win ratio):     ", formatC(x$log_win_ratio, digits = 4, format = "f"), "\n", sep = "")
  cat("Bootstrap SE (log): ", if (is.finite(x$se_log_win_ratio)) formatC(x$se_log_win_ratio, digits = 4, format = "f") else "NA", "\n", sep = "")
  cat(sprintf("%.1f%% CI for WR:    ", 100*x$conf_level),
      if (all(is.finite(x$ci_win_ratio))) paste0(formatC(x$ci_win_ratio[1], digits=4, format="f"), " - ", formatC(x$ci_win_ratio[2], digits=4, format="f")) else "NA",
      "\n", sep = "")
  cat("Bootstrap usable:   ", x$bootstrap$n_success, " / ", x$bootstrap$n_attempted,
      " (", formatC(100*x$bootstrap$success_fraction, digits=1, format="f"), "%)\n", sep = "")
  invisible(x)
}
