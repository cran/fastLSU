#' @title Fast Linear Step Up Procedure
#'
#' @description Fast Linear Step Up procedure of Benjamini–Hochberg FDR method

#' @param pvalues A vector of p values.
#' @param alpha The desired significance level, the default is 0.05.
#' @param algorithm Value of 1 or 2, 1 indicating algorithm 1 in the reference paper for a single chunk without sorting the p values; 2 indicating algorithm 2  in the reference paper for two or more chunks of p values.
#' @param mtg Global number of p values in the problem (the sum of all p values in all chunks), should be NA if algorithm is 1.
#' @param mtc Number of p values in current chunk (should be less than mtg), should be NA if algorithm is 1.
#' @export
#' @seealso \code{\link[stats]{p.adjust}}.
#' @references The "fastLSU" package is created based on Vered Madar and Sandra Batista’s work.  For more details:
#' @references Madar V, Batista S. FastLSU: a more practical approach for the Benjamini-Hochberg FDR controlling procedure for huge-scale testing problems. Bioinformatics 2016; 32:1716-23.
#' @return The result is a vector of significnat p values at the significance level of alpha.
#'
#' @examples
#' set.seed(111)
#' #simulate p-values for example;
#' pvals.sim = runif(30000)
#' Bsig = 1-rbinom(30000,1, .02)
#' Bsig[Bsig==0] = .0001
#' pvals.sim = pvals.sim*Bsig
#'
#' # Example of a single chunk
#' results.allchunks = fastLSU(pvalues=pvals.sim,alpha=0.1, algorithm = 1)
#'
#' # Example of 2 chunks of p-values (1st of 10,000, 2nd of 20,000)
#' pv.chunk1 = pvals.sim[1:10000]
#' pv.chunk2 = pvals.sim[10001:30000]
#' # Step 1:
#' results.chunk1 = fastLSU(pvalues=pv.chunk1,alpha=0.1, algorithm = 2,
#'  mtc=length(pv.chunk1),mtg=length(pvals.sim))
#'
#' # Step 2:
#' results.chunk2 = fastLSU(pvalues=pv.chunk2,alpha=0.1, algorithm = 2,
#' mtc=length(pv.chunk2),mtg=length(pvals.sim))
#'
#' # Step 3: # get final candidancy
#' cand.pvals = c(results.chunk1,results.chunk2)
#' # keep mtc equal to mtg, the fastLSU function already considers length(cand.pvals)

#' results.final = fastLSU(pvalues=cand.pvals,alpha=0.1, algorithm = 2,
#'  mtc=length(pvals.sim),mtg=length(pvals.sim))


fastLSU <- function(pvalues, alpha = 0.05, algorithm = 1, mtg = NA, mtc =NA) {
  pvalues <- pvalues[!is.na(pvalues)]
  if (min(pvalues) < 0 || max(pvalues) > 1) {
    stop("p values not in valid range [0, 1].")
  }
  if (alpha > 1 || alpha < 0) stop("alpha must be a single number between 0 and 1", call. = FALSE)
  if (!algorithm %in% c(1,2)) stop("algorithm must be a single number of 1 or 2", call. = FALSE)

  m = length(pvalues)

  if (algorithm == 1) {
    cat("FastLSU for a single chunk of p values.\n")
    kmax.r = m
    kmax.rb = kmax.r + 1
    stp = 0
    repeat {
      stp = stp + 1
      kmax = which(pvalues < alpha * kmax.r/m)
      pvalues = pvalues[kmax]
      kmax.r = length(kmax)
      if (kmax.r == kmax.rb) {break} else{
        kmax.rb = kmax.r
        cat(paste("Step ",stp,":",length(kmax)," p-values were selected for candidance of significant.\n"))
      }
    }
    return(pvalues)
  }

  if (algorithm == 2) {

    if (is.na(mtg) || is.na(mtc)) stop("mtg or mtc should not be NA.", call. = FALSE)

    if (m > mtc) stop("The input chunk size is ",m,", bigger than the original chunk size of ",mtc,call. = FALSE)
    if (mtc > mtg) stop("The chunk size is ",mtc,", bigger than the global number of p values of ",mtc,call. = FALSE)


    cat("FastLSU for two or more chunks of p values.\n")
    cat(m," p values are tested and the input chunk size is ",mtc,"\n")

    kmax.r = m
    kmax.rb = kmax.r + 1
    stp = 0
    repeat {
      stp = stp + 1
      kmax = which(pvalues < (kmax.r + mtg - mtc) * alpha / mtg)
      pvalues = pvalues[kmax]
      kmax.r = length(kmax)
      if (kmax.r == kmax.rb) {break}  else {
        kmax.rb = kmax.r
        cat(paste("Step ",stp,":",length(kmax)," p-values were selected for candidance of significant.\n"))
      }
    }
    return(pvalues)
  }
}
