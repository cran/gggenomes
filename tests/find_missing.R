library(tools)
pkg_of_interest <- "gggenomes"
db <- Rd_db(pkg_of_interest)

# Extract relevant fields / summary of fields into a database. This is from the example code on the ?Rd_db page. There can only be one `name` or `value` tag, but several `alias` tags

name <- lapply(db, tools:::.Rd_get_metadata, "name")
alias <- lapply(db, tools:::.Rd_get_metadata, "alias")
value <- lapply(db, tools:::.Rd_get_metadata, "value")

n_aliases <- lengths(alias)
df <- data.frame(
        file_name = rep(names(db), n_aliases),
        name = rep(unlist(name, use.names = FALSE), n_aliases),
        alias = unlist(alias, use.names = FALSE),
        has_value = rep(lengths(value) > 0, n_aliases)
)

# Create subsets of the database, and find the aliases that have no values. This is trying to allow for the possibility that an alias occurs in more than one help file (is this allowed?)

alias_with_value <- subset(df, has_value)
alias_without_value <- subset(df, !has_value)
no_value <- subset(alias_without_value, !alias %in% alias_with_value$alias)

# Find all the exports in the package, and subset the help pages to just those.

exports <- getNamespaceExports(pkg_of_interest)
subset(no_value, alias %in% exports)
