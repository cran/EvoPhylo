#' Post-process posterior tree samples from Bayesian phylogenetic analyses
#'
#' @description
#' An R implementation of the most commonly used functions of the Perl script
#' \href{https://github.com/nylander/Burntrees}{burntrees.pl} (J. Nylander).
#' \code{burntrees()} reads one or more posterior tree files (e.g., MrBayes
#' \code{.t} files from independent runs), removes burn-in, thins the samples,
#' optionally converts MrBayes clock trees from substitutions per site to time
#' units, strips comments, branch lengths and/or support values, and combines
#' everything into a single Nexus or Newick tree file and/or R object.
#'
#' @param path A character vector of tree file names, or the path to a single
#'   folder. If a folder is supplied, all files ending in \code{.t} (MrBayes)
#'   or \code{.trees} (BEAST2) in that folder are used, in alphabetical order.
#' @param burnin Either the number (\code{>= 1}) or the proportion
#'   (\code{< 1}) of trees to drop from the beginning of \emph{each} file, as
#'   in \code{\link{combine_log}}. Default is \code{0} (no burn-in).
#'   Equivalent to \code{burntrees.pl --burnin} / \code{--pburnin}.
#' @param thin A positive integer. After burn-in, every \code{thin}-th tree of
#'   each file is retained, starting with the first post-burn-in tree. Default
#'   is \code{1} (no thinning). Equivalent to \code{burntrees.pl --jump}.
#' @param myr \code{logical}; whether to transform the branch lengths of
#'   MrBayes clock trees from substitutions per site to time units by dividing
#'   each branch length by the clock rate stored in the
#'   \code{[&clockrate=...]} comment of each tree. Equivalent to
#'   \code{burntrees.pl --myr}. Default is \code{FALSE}.
#' @param rm.comments \code{logical}; whether to remove all comments
#'   (text within square brackets, e.g. \code{[&R]}, \code{[&clockrate=...]},
#'   \code{[&B IgrBrlens ...]}) from the trees. Equivalent to
#'   \code{burntrees.pl --rmcomments}. Default is \code{FALSE}.
#' @param rm.brlens \code{logical}; whether to remove branch lengths.
#'   Equivalent to \code{burntrees.pl --rmbrlens}. Default is \code{FALSE}.
#' @param rm.support \code{logical}; whether to remove numerical support
#'   values (node labels such as posterior probabilities or bootstrap values).
#'   Equivalent to \code{burntrees.pl --rmsupport}. Default is \code{FALSE}.
#' @param topology.only \code{logical}; shortcut that sets \code{rm.comments},
#'   \code{rm.brlens} and \code{rm.support} to \code{TRUE}, keeping only the
#'   taxa and tree topology. Equivalent to \code{burntrees.pl -R}. Default is
#'   \code{FALSE}.
#' @param trees.only \code{logical}; whether to output the tree lines only,
#'   without the Nexus header, translation table, and trailing
#'   \code{"end;"}. Equivalent to \code{burntrees.pl --treesonly}. Ignored
#'   when \code{format = "newick"}, which always contains trees only.
#'   Default is \code{FALSE}.
#' @param format The output format, one of \code{"nexus"} (default; the
#'   original Nexus format, keeping the translation table and taxon numbers if
#'   present), \code{"altnexus"} (Nexus format with taxon names in the trees
#'   and no translation table), or \code{"newick"} (one Newick tree per line,
#'   with taxon names; \code{"phylip"} is accepted as a synonym).
#'   Equivalent to \code{burntrees.pl --format}.
#' @param file An optional path to the output file. If \code{NULL} (default),
#'   no file is written.
#' @param output What to return in R: \code{"text"} (default) returns the
#'   processed trees as a character vector of lines in the chosen
#'   \code{format}; \code{"multiPhylo"} returns an \code{ape} \code{multiPhylo}
#'   object (comments are always dropped from these, see Details).
#'
#' @details
#' \code{burntrees()} processes tree files as text, line by line, in the same
#' way as \code{burntrees.pl}. This makes it fast for large posterior samples
#' and preserves all tree annotations (e.g., relaxed clock rates stored as
#' comments) unless they are explicitly removed. Each tree must be written on a
#' single line, as is the case for MrBayes, BEAST/BEAST2, and RevBayes output.
#'
#' Operations are applied in the following order: burn-in removal and thinning
#' (both applied separately to each file), conversion to time units
#' (\code{myr}), translation of taxon numbers into names (for
#' \code{format = "altnexus"} or \code{"newick"}), removal of support values,
#' branch lengths, and comments, and finally combination of all files.
#'
#' \strong{Combining runs.} When several files are supplied, the trees
#' retained from each file are combined, in the order of \code{path}, into a
#' single tree set. With \code{format = "nexus"}, the header and translation
#' table of the first file are used. If the translation tables differ among
#' files (but contain the same taxa), the taxon numbers of all trees are
#' re-mapped to the translation table of the first file.
#'
#' \strong{Clock trees in time units.} In MrBayes clock trees, branch lengths
#' are in expected substitutions per site, i.e., time multiplied by the base
#' clock rate. Each tree in a MrBayes clock \code{.t} file stores the sampled
#' clock rate in a \code{[&clockrate=...]} comment. With \code{myr = TRUE},
#' each branch length is divided by that tree's clock rate, so branch lengths
#' are returned in the same time units used for the calibrations (typically
#' millions of years). Values are printed with six decimal places, as in
#' \code{burntrees.pl}. Other annotations (e.g., \code{[&B IgrBrlens ...]})
#' are left unchanged.
#'
#' \strong{multiPhylo output.} With \code{output = "multiPhylo"}, the
#' processed trees are read with \code{\link[ape]{read.tree}} after
#' translating taxon numbers into names and removing comments. Tree names are
#' taken from the tree file (e.g., \code{"gen.10000"}), prefixed with the file
#' name when several files are combined. If \code{file} is also supplied, the
#' file is written in the chosen \code{format}, with comments kept unless
#' \code{rm.comments = TRUE}.
#'
#' Unlike \code{burntrees.pl}, posterior probabilities stored as comments in
#' MrBayes consensus trees (\code{.con.tre}) are not converted into support
#' values when comments are removed; use \code{treeio::read.mrbayes()} for
#' those files.
#'
#' @return
#' If \code{output = "text"}, a character vector with one element per line of
#' the processed tree file. If \code{output = "multiPhylo"}, an object of
#' class \code{multiPhylo}. If \code{file} is supplied, the output is written
#' to \code{file} and the object is returned invisibly.
#'
#' @author J. A. A. Nylander (original Perl script \code{burntrees.pl});
#'   R implementation by T. R. Simões.
#'
#' @seealso \code{\link{combine_log}} for the equivalent processing of
#'   parameter log (\code{.p}) files.
#'
#' @examples
#' # Two runs of a MrBayes tip-dated relaxed clock (IGR) analysis
#' run1 <- system.file("extdata", "ex_mb_clock.run1.t", package = "EvoPhylo")
#' run2 <- system.file("extdata", "ex_mb_clock.run2.t", package = "EvoPhylo")
#'
#' # Remove 25% burn-in from each run, keep every 5th tree, and combine
#' # the two runs into a single MrBayes-style Nexus tree set
#' trees <- burntrees(c(run1, run2), burnin = 0.25, thin = 5)
#' writeLines(substr(trees[1:20], 1, 70))
#'
#' # Convert branch lengths to time units (Myr) and drop all comments
#' trees_myr <- burntrees(c(run1, run2), burnin = 0.25, thin = 5,
#'                        myr = TRUE, rm.comments = TRUE)
#'
#' # Same, returned as a multiPhylo object with taxon names
#' trees_phylo <- burntrees(c(run1, run2), burnin = 0.25, thin = 5,
#'                          myr = TRUE, output = "multiPhylo")
#' trees_phylo
#'
#' # Topologies only, as Newick strings with taxon names
#' burntrees(run1, burnin = 0.25, thin = 10, topology.only = TRUE,
#'           format = "newick")
#'
#' \dontrun{
#' # Write the combined, time-scaled trees to a Nexus file
#' burntrees(c(run1, run2), burnin = 0.25, myr = TRUE,
#'           file = "combined_myr.t")
#' }
#'
#' @export
burntrees <- function(path, burnin = 0, thin = 1, myr = FALSE,
                      rm.comments = FALSE, rm.brlens = FALSE,
                      rm.support = FALSE, topology.only = FALSE,
                      trees.only = FALSE,
                      format = c("nexus", "altnexus", "newick"),
                      file = NULL, output = c("text", "multiPhylo")) {

  ## Check arguments ---------------------------------------------------------
  if (!is.character(path) || length(path) == 0) {
    files <- NULL
  } else if (all(utils::file_test("-f", path))) {
    files <- path
  } else if (length(path) == 1 && utils::file_test("-d", path)) {
    files <- list.files(path, pattern = "\\.(t|trees)$", full.names = TRUE)
  } else {
    files <- NULL
  }
  if (length(files) == 0) {
    stop("'path' must be a character vector containing the names of tree files or of a folder containing such files (.t or .trees).",
         call. = FALSE)
  }

  if (length(burnin) != 1 || !is.numeric(burnin) || is.na(burnin) || burnin < 0) {
    stop("'burnin' must be a single number corresponding to the number or proportion of trees to drop from the beginning of each tree file.",
         call. = FALSE)
  }
  if (length(thin) != 1 || !is.numeric(thin) || is.na(thin) || thin < 1 ||
      thin != round(thin)) {
    stop("'thin' must be a single positive integer; every 'thin'-th tree is retained.",
         call. = FALSE)
  }

  for (a in c("myr", "rm.comments", "rm.brlens", "rm.support",
              "topology.only", "trees.only")) {
    v <- get(a)
    if (length(v) != 1 || !is.logical(v) || is.na(v)) {
      stop(sprintf("'%s' must be TRUE or FALSE.", a), call. = FALSE)
    }
  }

  if (!is.character(format)) {
    stop("'format' must be one of \"nexus\", \"altnexus\", or \"newick\".", call. = FALSE)
  }
  format <- tolower(format[1])
  if (format == "phylip") format <- "newick"
  format <- match.arg(format, c("nexus", "altnexus", "newick"))
  output <- match.arg(output)

  if (!is.null(file) && (!is.character(file) || length(file) != 1)) {
    stop("'file' must be NULL or a single file name.", call. = FALSE)
  }

  if (topology.only) {
    rm.comments <- rm.brlens <- rm.support <- TRUE
  }
  if (myr && rm.brlens) {
    warning("Branch lengths are removed, so 'myr = TRUE' has no effect.",
            call. = FALSE)
    myr <- FALSE
  }

  ## Read files, remove burn-in and thin ---------------------------------------
  L <- lapply(files, .bt_read_treefile)
  names(L) <- files

  L <- lapply(L, function(x) {
    n <- length(x$desc)
    if (burnin > 0) {
      if (burnin < 1) {
        b <- seq_len(round(burnin * n))
      } else {
        b <- seq_len(round(min(burnin, n)))
      }
      if (length(b) > 0) {
        x$desc <- x$desc[-b]
        x$names <- x$names[-b]
      }
    }
    if (length(x$desc) == 0) {
      stop(sprintf("No trees left in file '%s' after removing the burn-in.", x$file),
           call. = FALSE)
    }
    keep <- seq(1, length(x$desc), by = thin)
    x$desc <- x$desc[keep]
    x$names <- x$names[keep]
    x
  })

  ## Taxon numbers or taxon names in the output? ------------------------------
  ref <- L[[1]]$translate
  keyed <- format == "nexus" && !is.null(ref)

  L <- lapply(L, function(x) {
    ## Convert to time units
    if (myr) x$desc <- .bt_brlen2time(x$desc, x$file)

    ## Taxon translation
    tr <- x$translate
    if (keyed) {
      if (!is.null(tr) && !identical(tr, ref)) {
        inv <- stats::setNames(names(ref), ref)
        map <- stats::setNames(inv[tr], names(tr))
        if (anyNA(map)) {
          stop(sprintf("Taxa %s in file '%s' are not found in the translation table of the first file.",
                       paste(tr[is.na(map)], collapse = ", "), x$file), call. = FALSE)
        }
        x$desc <- .bt_map_tips(x$desc, map)
      } else if (is.null(tr)) {
        x$desc <- .bt_map_tips(x$desc, stats::setNames(names(ref), ref))
      }
    } else if (!is.null(tr)) {
      x$desc <- .bt_map_tips(x$desc, tr)
    }

    ## Remove support values, branch lengths, and comments
    if (rm.support) x$desc <- gsub(.bt_re$support, "", x$desc, perl = TRUE)
    if (rm.brlens) x$desc <- gsub(.bt_re$brlen, "", x$desc, perl = TRUE)
    if (rm.comments) {
      x$desc <- trimws(gsub(.bt_re$comment, "", x$desc, perl = TRUE))
      x$names <- trimws(gsub(.bt_re$comment, "", x$names, perl = TRUE))
    }
    x
  })

  ## Combine files --------------------------------------------------------------
  desc <- unlist(lapply(L, `[[`, "desc"), use.names = FALSE)
  tree_names <- unlist(lapply(L, `[[`, "names"), use.names = FALSE)

  ## Build output text ----------------------------------------------------------
  text <- NULL
  if (output == "text" || !is.null(file)) {
    if (format == "newick") {
      text <- desc
    } else {
      tree_lines <- paste0("   tree ", tree_names, " = ", desc)
      if (trees.only) {
        text <- tree_lines
      } else {
        header <- L[[1]]$header
        if (length(header) == 0) {
          header <- c("#NEXUS", "begin trees;")
        } else if (!keyed) {
          header <- .bt_drop_translate(header)
        }
        text <- c(header, tree_lines, "end;")
      }
    }
  }

  if (!is.null(file)) {
    writeLines(text, file)
  }

  ## Return -----------------------------------------------------------------------
  if (output == "text") {
    if (is.null(file)) return(text) else return(invisible(text))
  }

  phy_desc <- if (keyed) .bt_map_tips(desc, ref) else desc
  phy_desc <- trimws(gsub(.bt_re$comment, "", phy_desc, perl = TRUE))
  phy_names <- trimws(gsub(.bt_re$comment, "", tree_names, perl = TRUE))
  if (length(files) > 1) {
    run_id <- rep(tools::file_path_sans_ext(basename(files)),
                  vapply(L, function(x) length(x$desc), integer(1)))
    phy_names <- paste(run_id, phy_names, sep = "_")
  }
  trees <- ape::read.tree(text = phy_desc, tree.names = phy_names)
  if (inherits(trees, "phylo")) {
    trees <- structure(list(trees), class = "multiPhylo")
    names(trees) <- phy_names
  }

  if (is.null(file)) trees else invisible(trees)
}


## Internal helpers for burntrees() -------------------------------------------------

## Regular expressions. Patterns starting with .bt_skip ignore anything inside
## square-bracket comments, e.g. "[&prob_range={1.0,1.0}]" or
## "[&B IgrBrlens 5.39e-01]".
.bt_num <- "[-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?"
.bt_skip <- "\\[[^\\]]*\\](*SKIP)(*FAIL)|"
.bt_re <- list(
  comment   = "\\[[^\\]]*\\]",
  brlen     = paste0(.bt_skip, ":\\s*", .bt_num),
  support   = paste0(.bt_skip, "(?<=\\))\\s*", .bt_num),
  tip       = paste0(.bt_skip, "[(,]\\s*\\K(?:'(?:[^']|'')*'|[^\\s(),:;\\[\\]']+)"),
  clockrate = paste0("\\[&clockrate\\s*=\\s*(", .bt_num, ")\\s*\\]"),
  tree      = "^\\s*tree\\s+\\*?\\s*((?:[^=\\[]|\\[[^\\]]*\\])+?)\\s*=\\s*(.*)$"
)

## Read a Nexus or Newick tree file. Returns a list with the file name, the
## header (Nexus lines before the first tree), the translation table (named
## character vector: names = numbers/keys, values = taxon names, or NULL),
## the tree names, and the tree descriptions (everything after "tree name =").
.bt_read_treefile <- function(f) {
  lines <- readLines(f, warn = FALSE)
  lines <- sub("\r$", "", lines)
  first <- lines[nzchar(trimws(lines))][1]

  if (!is.na(first) && grepl("^\\s*#nexus", first, ignore.case = TRUE)) {
    begin <- grep("^\\s*begin\\s+trees\\s*;", lines, ignore.case = TRUE)
    tree_idx <- grep("^\\s*tree\\s", lines, ignore.case = TRUE)
    if (length(begin) > 0) tree_idx <- tree_idx[tree_idx > begin[1]]
    if (length(tree_idx) == 0) {
      stop(sprintf("No trees found in file '%s'.", f), call. = FALSE)
    }
    header <- lines[seq_len(tree_idx[1] - 1)]
    translate <- .bt_translate(header, f)

    m <- regmatches(lines[tree_idx],
                    regexec(.bt_re$tree, lines[tree_idx], ignore.case = TRUE, perl = TRUE))
    ok <- lengths(m) == 3
    if (!all(ok)) {
      stop(sprintf("Could not read tree on line %s of file '%s'.",
                   tree_idx[!ok][1], f), call. = FALSE)
    }
    tree_names <- trimws(vapply(m, `[`, character(1), 2))
    desc <- trimws(vapply(m, `[`, character(1), 3))
  } else {
    tree_idx <- grep("\\(.*;\\s*$", lines)
    if (length(tree_idx) == 0) {
      stop(sprintf("File '%s' does not seem to be a Nexus or Newick tree file.", f),
           call. = FALSE)
    }
    header <- character(0)
    translate <- NULL
    desc <- trimws(lines[tree_idx])
    tree_names <- paste0("tree_", seq_along(desc))
  }

  if (!all(grepl(";$", desc))) {
    stop(sprintf("Trees in file '%s' must be written on a single line each.", f),
         call. = FALSE)
  }

  list(file = f, header = header, translate = translate,
       names = tree_names, desc = desc)
}

## Parse the translation table from the Nexus lines preceding the first tree
.bt_translate <- function(header, f) {
  i <- grep("^\\s*translate\\b", header, ignore.case = TRUE, perl = TRUE)
  if (length(i) == 0) return(NULL)
  txt <- paste(header[i[1]:length(header)], collapse = " ")
  txt <- sub("^\\s*translate\\b", "", txt, ignore.case = TRUE, perl = TRUE)
  txt <- sub(";.*$", "", txt)
  entries <- strsplit(txt, ",(?=(?:[^']*'[^']*')*[^']*$)", perl = TRUE)[[1]]
  entries <- trimws(entries)
  entries <- entries[nzchar(entries)]
  key <- sub("^(\\S+)\\s+.*$", "\\1", entries)
  lab <- trimws(sub("^\\S+\\s+", "", entries))
  if (anyDuplicated(key) || any(key == lab & !grepl("\\s", entries))) {
    stop(sprintf("Could not read the translation table of file '%s'.", f),
         call. = FALSE)
  }
  stats::setNames(lab, key)
}

## Remove the translate block from a Nexus header
.bt_drop_translate <- function(header) {
  i <- grep("^\\s*translate\\b", header, ignore.case = TRUE, perl = TRUE)
  if (length(i) == 0) return(header)
  j <- grep(";", header)
  j <- j[j >= i[1]][1]
  header[-(i[1]:j)]
}

## Vectorised "gsub with a function": all matches of 'pattern' in 'x' are
## passed at once to fun(tokens, index_of_string), which returns replacements
.bt_gsub_fun <- function(x, pattern, fun) {
  m <- gregexpr(pattern, x, perl = TRUE)
  pieces <- regmatches(x, m, invert = NA)    # gap, match, gap, match, ..., gap
  n <- lengths(pieces)
  flat <- unlist(pieces, use.names = FALSE)
  is_tok <- sequence(n) %% 2L == 0L
  if (!any(is_tok)) return(x)
  id <- rep.int(seq_along(x), n)
  flat[is_tok] <- fun(flat[is_tok], id[is_tok])
  ends <- cumsum(n)
  starts <- ends - n + 1L
  vapply(seq_along(x), function(i) paste(flat[starts[i]:ends[i]], collapse = ""),
         character(1))
}

## Replace tip tokens (taxon numbers or names) using a named character vector
.bt_map_tips <- function(desc, map) {
  .bt_gsub_fun(desc, .bt_re$tip, function(tok, id) {
    new <- map[tok]
    new[is.na(new)] <- tok[is.na(new)]
    unname(new)
  })
}

## Divide all branch lengths of each tree by its [&clockrate=...] value
.bt_brlen2time <- function(desc, f) {
  cr <- regmatches(desc, regexec(.bt_re$clockrate, desc, ignore.case = TRUE, perl = TRUE))
  rate <- suppressWarnings(as.numeric(vapply(cr, function(x) {
    if (length(x) == 2) x[2] else NA_character_
  }, character(1))))
  if (anyNA(rate)) {
    stop(sprintf("No [&clockrate=...] comment found in tree %s of file '%s'. 'myr = TRUE' only applies to MrBayes clock trees.",
                 which(is.na(rate))[1], f), call. = FALSE)
  }
  if (any(rate <= 0)) {
    stop(sprintf("Clock rate <= 0 found in tree %s of file '%s'.",
                 which(rate <= 0)[1], f), call. = FALSE)
  }
  .bt_gsub_fun(desc, .bt_re$brlen, function(tok, id) {
    sprintf(":%f", as.numeric(substring(tok, 2)) / rate[id])
  })
}
