#' Statistical tests to compare word clouds
#'
#' \code{\link{wordcloudstats}} is the main user interface to conduct a quantitative analysis of
#' the significance of the observed differences between the two word clouds.
#'
#' @name wordcloudanalysis-package
#' @docType package
NULL

#' @export
#' @title Quantitative statistical analysis of differences between two word clouds.
#' @description Main user interface function for quantitative statistical analysis between two
#' word clouds. If \code{avoid.words} is specified, then words that ought to be ignored are
#' excluded from the analysis. The common terms between \code{l1} and \code{l2} are
#' considered for statistical analysis. The analysis computes the significance of the
#' observed differences in term occurrence in the \code{l1} and \code{l2}.
#'
#' @param l1 A character vector representing terms corresponding to one word cloud
#' @param l2 A character vector representing terms corresponding to the other word cloud
#' @param names A character vector of length 2 with the two names for the provided lists. The
#'              first element corresponds to \code{l1} and the second element corresponds to \code{l2}.
#' @param avoid.words A character vector of words to ignore when analyzing the two word clouds.
#'                    Defaults to \code{NULL}.
#'
#' @return A list of three elements. \code{outputs} contains a matrix of the significance and
#' counts of each of the common terms between the lists. \code{counts} contains a table of
#' the counts of each term. \code{frequency} contains a table of the frequencies of each term
#' in each list.
#'
#' @examples
#' group1 <- c('head','toe','hand',rep('knee',4))
#' group2 <- c(rep('toe',3),'hand',rep('head',2))
#'
#' stats <- wordcloudstats(group1, group2, names=c('Group 1','Group 2'))
#'
#' head(stats$outputs)
#'
#' @seealso \code{\link[stats]{prop.test}}
wordcloudstats <- function(l1, l2, names, avoid.words=NULL) {

  # Remove words that the user wants to avoid
  if (! is.null(avoid.words)) {
    l1 <- l1[!l1 %in% avoid.words]
    l2 <- l2[!l2 %in% avoid.words]
  }

  # Keep only the terms that appear at least once in each list
  common <- intersect(unique(l1),unique(l2))
  l1 <- l1[l1%in%common]
  l2 <- l2[l2%in%common]

  # Basic frequency table
  l1.freq <- unlist(table(l1))/length(l1)
  l2.freq <- unlist(table(l2))/length(l2)

  freq <- matrix(NA,nrow=length(l1.freq),ncol=2)
  rownames(freq) <- names(l1.freq)
  freq[,1] <- l1.freq
  freq[names(l2.freq),2] <- l2.freq
  colnames(freq) <- names

  # Basic count table
  l1.ct <- unlist(table(l1))
  l2.ct <- unlist(table(l2))

  count <- matrix(NA,nrow=length(l1.ct),ncol=2)
  rownames(count) <- names(l1.ct)
  count[,1] <- l1.ct
  count[names(l2.ct),2] <- l2.ct
  colnames(count) <- names

  pvals <- rep(1,nrow(count))
  names(pvals) <- rownames(count)
  g.in <- rep(1,nrow(count))
  names(g.in) <- rownames(count)

  l1.in <- rep(1,nrow(count))
  names(l1.in) <- rownames(count)
  l1.out <- rep(1,nrow(count))
  names(l1.out) <- rownames(count)

  l2.in <- rep(1,nrow(count))
  names(l2.in) <- rownames(count)
  l2.out <- rep(1,nrow(count))
  names(l2.out) <- rownames(count)

  # Run proportions test
  for (term in rownames(count)) {
    l1.yes <- count[term,1]
    l2.yes <- count[term,2]
    l1.no <- length(l1) - l1.yes
    l2.no <- length(l2) - l2.yes

    l1.in[term] <- l1.yes
    l1.out[term] <- l1.no
    l2.in[term] <- l2.yes
    l2.out[term] <- l2.no

    mat <- matrix(NA, nrow=2, ncol=2)
    mat[1,1] <- l1.yes
    mat[1,2] <- l1.no
    mat[2,1] <- l2.yes
    mat[2,2] <- l2.no

    ptest <- prop.test(mat)
    pvals[term] <- ptest$p.value
    maxval<- max(as.numeric(unlist(ptest$estimate)))
    g.in.ind <- grep(maxval,as.numeric(unlist(ptest$estimate)))
    g.in[term] <- names[g.in.ind]
  }

  p.adj <- p.adjust(unname(pvals),method='fdr')
  mat.val <- cbind(l1.in,l1.out,l2.in,l2.out,pvals,p.adj,g.in)
  present.names <- paste('Number of occurences in',names)
  absent.names <- paste('Number of unmatched terms in',names)
  mat.names <- c(present.names[1],absent.names[1], present.names[2],absent.names[2],'P-Value','BH Value','Greater frequency in')
  colnames(mat.val) <- mat.names

  list(outputs=mat.val,counts=count,frequency=freq)
}
