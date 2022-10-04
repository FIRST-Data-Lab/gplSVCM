ZZ.func <- function(Tr.set, ZZ.tilde, p, index.tr) {
  idx <- unlist(lapply(Tr.set, function(x) if (index.tr %in% x) which(x == index.tr) else 0))
  all.ZZ <- matrix(0, p, p)
  for (iter in 1:length(idx)) {
    start <- (idx[iter] - 1) * p + 1
    end <- (idx[iter]) * p
    if (idx[iter] != 0) {
      temp <- ZZ.tilde[[iter]][start:end, start:end]
      all.ZZ <- all.ZZ + temp
    }
  }
  sum(idx != 0) * all.ZZ
}
