# Helper functions for dog_mutation_project.R
#
# These functions replace the machine-specific `common.R` file used during
# exploratory development. They use only base R and lme4/lmerTest objects.

# Qualitative palettes used throughout the manuscript figures.
cs_gg <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#FFFF33", "#A65628", "#F781BF", "#999999"
)
cs25 <- grDevices::hcl.colors(25, palette = "Dynamic")

# Split an ordered vector into approximately equal-sized groups.
chunk <- function(x, n) {
    split(x, factor(sort(rank(x) %% n)))
}

# Cosine similarity between two numeric vectors.
cos_sim <- function(x, y) {
    denominator <- sqrt(sum(x^2)) * sqrt(sum(y^2))
    if (denominator == 0) {
        return(NA_real_)
    }
    sum(x * y) / denominator
}

# Draw a fixed-effect regression line and optional 95% confidence band from a
# fitted mixed-effects model. The `species` argument is retained for backward
# compatibility with the original analysis calls.
fun_addRegressionLine_lmer <- function(model, variable, x_min, x_max,
                                       species = NULL, se = TRUE, ...) {
    x <- seq(x_min, x_max, length.out = 200)
    new_data <- data.frame(x)
    names(new_data) <- variable

    fixed_formula <- lme4::nobars(stats::formula(model))
    design <- stats::model.matrix(stats::delete.response(stats::terms(fixed_formula)), new_data)
    coefficients <- lme4::fixef(model)
    design <- design[, names(coefficients), drop = FALSE]
    fitted <- drop(design %*% coefficients)

    dots <- list(...)
    line_colour <- if (!is.null(dots$col)) dots$col else graphics::par("fg")

    if (isTRUE(se)) {
        covariance <- as.matrix(stats::vcov(model))
        standard_error <- sqrt(rowSums((design %*% covariance) * design))
        graphics::polygon(
            c(x, rev(x)),
            c(fitted - 1.96 * standard_error, rev(fitted + 1.96 * standard_error)),
            col = grDevices::adjustcolor(line_colour, alpha.f = 0.20),
            border = NA
        )
    }

    dots$col <- line_colour
    do.call(graphics::lines, c(list(x = x, y = fitted), dots))
    invisible(data.frame(x = x, fitted = fitted))
}

# Boxplot wrapper used for signature-robustness panels.
fun_boxplot <- function(x, points = FALSE, cex = 1, ...) {
    graphics::boxplot(x, outline = FALSE, xaxt = "n", ...)
    if (isTRUE(points)) {
        for (i in seq_along(x)) {
            graphics::points(rep(i, length(x[[i]])), x[[i]], pch = 20, cex = cex)
        }
    }
    invisible(x)
}

# Add a log10 axis labelled on the original scale.
add_axis <- function(side, at = 1:7) {
    graphics::axis(side, at = at, labels = format(10^at, scientific = TRUE))
}

fun_snv_profile <- function(profile, main = "", ylab = "%") {
    values <- profile / sum(profile) * 100
    colours <- rep(cs_gg[1:6], each = 16)
    graphics::barplot(
        values, col = colours, border = "white", space = 0,
        names.arg = names(profile), las = 3, main = main, ylab = ylab
    )
    graphics::abline(v = (1:5) * 16)
    graphics::text(
        (1:6) * 16 - 8, max(values) * 19 / 20,
        labels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    )
    invisible(values)
}

fun_indel_profile <- function(profile, main = "", ylab = "n. sindel") {
    colours <- grDevices::hcl.colors(length(profile), palette = "Spectral")
    graphics::barplot(
        profile, col = colours, border = "white", space = 0,
        names.arg = names(profile), las = 3, main = main, ylab = ylab
    )
    invisible(profile)
}
