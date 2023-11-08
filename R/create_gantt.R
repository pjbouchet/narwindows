#' Generate Gantt chat
#'
#' @export
#' 
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' }

create_gantt <- function(project, spots = NULL, 
                         by_date = FALSE, 
                         exact_date = FALSE,
                         project_start_date = Sys.Date(),
                         project_end_date = Sys.Date() + 365,
                         colour_palette = wesanderson::wes_palette("Darjeeling1"),
                         font_family = "sans", 
                         mark_quarters = FALSE, 
                         mark_years = FALSE,
                         size_wp = 6, 
                         hide_wp = FALSE, 
                         size_activity = 4, 
                         size_text_relative = 1,
                         label_wrap = FALSE, 
                         month_number_label = TRUE, 
                         month_date_label = TRUE,
                         x_axis_position = "top", 
                         colour_stripe = "lightgray", 
                         alpha_wp = 1,
                         alpha_activity = 1, 
                         line_end = "round", 
                         month_breaks = 1,
                         show_vertical_lines = TRUE, 
                         axis_text_align = "right") {
  
  # Adapted from ganttrify::ganttrify
  
  if (length(unique(project$wp)) > length(as.character(wesanderson::wes_palette("Darjeeling1")))) {
    colour_palette <- rep(colour_palette, length(unique(project$wp)))[1:length(unique(project$wp))]
  }
  if (label_wrap != FALSE) {
    if (isTRUE(label_wrap)) {
      label_wrap <- 32
    }
    project$wp <- stringr::str_wrap(
      string = project$wp,
      width = label_wrap
    )
    project$activity <- stringr::str_wrap(
      string = project$activity,
      width = label_wrap
    )
    if (is.null(spots) == FALSE) {
      spots$activity <- stringr::str_wrap(
        string = spots$activity,
        width = label_wrap
      )
    }
  }
  
  if (by_date == FALSE) {
    df <- project |> dplyr::mutate(
      start_date = as.numeric(start_date),
      end_date = as.numeric(end_date)
    )
    
    if (sum(is.na(df$start_date)) == nrow(df)) {
      stop("Make sure that the input data are properly formatted and/or you have selected the right paramaters.")
    }
    start_yearmon <- zoo::as.yearmon(project_start_date) -
      (1 / 12)
    df_yearmon <- df  |> 
      dplyr::mutate(start_date_yearmon = start_yearmon +
                      (1 / 12) * start_date, end_date_yearmon = start_yearmon +
                      (1 / 12) * zoo::as.yearmon(end_date)) |>
      dplyr::transmute(
        wp = as.character(wp),
        activity = as.character(activity), start_date = zoo::as.Date(start_date_yearmon,
                                                                     frac = 0
        ), end_date = zoo::as.Date(end_date_yearmon,
                                   frac = 1
        )
      )
  } else {
    if (exact_date == TRUE) {
    } else {
      df_yearmon <- project |>
        dplyr::mutate(
          start_date_yearmon = zoo::as.yearmon(start_date),
          end_date_yearmon = zoo::as.yearmon(end_date)
        ) |>
        dplyr::transmute(
          wp = as.character(wp), activity = as.character(activity),
          start_date = zoo::as.Date(start_date_yearmon,
                                    frac = 0
          ), end_date = zoo::as.Date(end_date_yearmon,
                                     frac = 1
          )
        )
    }
  }
  if (exact_date == TRUE) {
    df <- project |> dplyr::mutate(
      start_date = as.Date(start_date),
      end_date = as.Date(end_date), wp = as.character(wp),
      activity = as.character(activity)
    )
    df_yearmon <- df |> dplyr::mutate(start_date = zoo::as.Date(zoo::as.yearmon(start_date),
                                                                frac = 0
    ), end_date = zoo::as.Date(zoo::as.yearmon(end_date),
                               frac = 1
    ))
  }
  
  # Phil == Changed
  # sequence_months <- seq.Date(
  #   from = min(df_yearmon[["start_date"]]),
  #   to = max(df_yearmon[["end_date"]]), by = "1 month"
  # )
  
  sequence_months <- seq.Date(
    from = lubridate::as_date(paste0(min(lubridate::year(project$start_date)), "-01-01")),
    to = lubridate::as_date(paste0(max(lubridate::year(project$end_date)), "-12-31")), by = "1 month"
  )
  
  if (length(sequence_months) %% 2 != 0) {
    sequence_months <- seq.Date(
      from = min(df_yearmon[["start_date"]]),
      to = max(df_yearmon[["end_date"]]) + 1, by = "1 month"
    )
  }
  date_range_matrix <- matrix(as.numeric(sequence_months),
                              ncol = 2, byrow = TRUE
  )
  date_range_df <- tibble::tibble(start = zoo::as.Date.numeric(date_range_matrix[
    ,
    1
  ]), end = zoo::as.Date.numeric(date_range_matrix[, 2]))
  
  
  # Phil == Changed
  date_breaks <- zoo::as.Date(zoo::as.yearmon(seq.Date(from = lubridate::as_date(paste0(min(lubridate::year(project$start_date)), "-01-01")),
                                                       to = lubridate::as_date(paste0(max(lubridate::year(project$end_date)), "-12-31")), 
                                                       by = paste(month_breaks, "month"))), frac = 0.5)
  
  # date_breaks <- zoo::as.Date(zoo::as.yearmon(seq.Date(from = min(df_yearmon[["start_date"]] +
  #   15), to = max(df_yearmon[["end_date"]] + 15), by = paste(
  #   month_breaks,
  #   "month"
  # ))), frac = 0.5)
  
  date_breaks_q <- seq.Date(from = lubridate::floor_date(
    x = min(df_yearmon[["start_date"]]),
    unit = "year"
  ), to = lubridate::ceiling_date(
    x = max(df_yearmon[["end_date"]]),
    unit = "year"
  ), by = "1 quarter")
  
  date_breaks_y <- seq.Date(from = lubridate::floor_date(
    x = min(df_yearmon[["start_date"]]),
    unit = "year"
  ), to = lubridate::ceiling_date(
    x = max(df_yearmon[["end_date"]]),
    unit = "year"
  ), by = "1 year")
  df_levels <- rev(df_yearmon |> dplyr::select(wp, activity) |>
                     t() |> as.character() |> unique())
  if (exact_date == TRUE) {
    df_yearmon_fct <- dplyr::bind_rows(
      activity = df, wp = df |>
        dplyr::group_by(wp) |> dplyr::summarise(
          activity = unique(wp),
          start_date = min(start_date), end_date = max(end_date)
        ),
      .id = "type"
    ) |>
      dplyr::mutate(activity = factor(
        x = activity,
        levels = df_levels
      )) |>
      dplyr::arrange(activity)
  } else {
    df_yearmon_fct <- dplyr::bind_rows(
      activity = df_yearmon,
      wp = df_yearmon |> dplyr::group_by(wp) |> dplyr::summarise(
        activity = unique(wp),
        start_date = min(start_date), end_date = max(end_date)
      ),
      .id = "type"
    ) |>
      dplyr::mutate(activity = factor(
        x = activity,
        levels = df_levels
      )) |>
      dplyr::arrange(activity)
  }
  if (hide_wp == TRUE) {
    df_yearmon_fct <- df_yearmon_fct |> dplyr::filter(type !=
                                                        "wp")
  }
  gg_gantt <- ggplot2::ggplot(data = df_yearmon_fct, mapping = ggplot2::aes(
    x = start_date,
    y = activity, xend = end_date, yend = activity, colour = wp
  )) +
    ggplot2::geom_rect(
      data = date_range_df, ggplot2::aes(
        xmin = start,
        xmax = end, ymin = -Inf, ymax = Inf
      ), inherit.aes = FALSE,
      alpha = 0.4, fill = colour_stripe
    )
  if (mark_quarters == TRUE) {
    gg_gantt <- gg_gantt + ggplot2::geom_vline(
      xintercept = date_breaks_q,
      colour = "gray50"
    )
  }
  if (mark_years == TRUE) {
    gg_gantt <- gg_gantt + ggplot2::geom_vline(
      xintercept = date_breaks_y,
      colour = "gray50"
    )
  }
  df_yearmon_fct$wp_alpha <- 0
  df_yearmon_fct$activity_alpha <- 0
  df_yearmon_fct <- df_yearmon_fct |> dplyr::mutate(activity_alpha = ifelse(type ==
                                                                              "activity", alpha_activity, 0))
  df_yearmon_fct <- df_yearmon_fct |> dplyr::mutate(wp_alpha = ifelse(type ==
                                                                        "wp", alpha_wp, 0))
  gg_gantt <- gg_gantt + ggplot2::geom_segment(
    data = df_yearmon_fct,
    lineend = line_end, size = size_activity, alpha = df_yearmon_fct$activity_alpha
  ) +
    ggplot2::geom_segment(
      data = df_yearmon_fct, lineend = line_end,
      size = size_wp, alpha = df_yearmon_fct$wp_alpha
    )
  if (month_number_label == TRUE & month_date_label == TRUE) {
    gg_gantt <- gg_gantt + ggplot2::scale_x_date(
      name = "",
      breaks = date_breaks, date_labels = "%b\n", minor_breaks = NULL,
      sec.axis = ggplot2::dup_axis(labels = paste0(
        "M",
        seq_along(date_breaks) * month_breaks - (month_breaks -
                                                   1)
      ))
    )
  } else if (month_number_label == FALSE & month_date_label ==
             TRUE) {
    gg_gantt <- gg_gantt + ggplot2::scale_x_date(
      name = "",
      breaks = date_breaks, date_labels = "%b\n", minor_breaks = NULL,
      position = x_axis_position
    )
  } else if (month_number_label == TRUE & month_date_label ==
             FALSE) {
    gg_gantt <- gg_gantt + ggplot2::scale_x_date(
      name = "",
      breaks = date_breaks, date_labels = paste0("M", seq_along(date_breaks)),
      minor_breaks = NULL, position = x_axis_position
    )
  } else if (month_number_label == FALSE & month_date_label ==
             FALSE) {
    gg_gantt <- gg_gantt + ggplot2::scale_x_date(name = "")
  }
  if (axis_text_align == "right") {
    axis_text_align_n <- 1
  } else if (axis_text_align == "centre" | axis_text_align ==
             "center") {
    axis_text_align_n <- 0.5
  } else if (axis_text_align == "left") {
    axis_text_align_n <- 0
  } else {
    axis_text_align_n <- 1
  }
  gg_gantt <- suppressWarnings(gg_gantt + ggplot2::scale_y_discrete("") +
                                 ggplot2::theme_minimal() + ggplot2::scale_colour_manual(values = colour_palette) +
                                 ggplot2::theme(
                                   text = ggplot2::element_text(family = font_family),
                                   axis.text.y.left = ggplot2::element_text(
                                     face = ifelse(test = df_yearmon_fct |>
                                                     dplyr::distinct(activity, wp, type) |> dplyr::pull(type) ==
                                                     "wp", yes = "bold", no = "plain"), size = ggplot2::rel(size_text_relative),
                                     hjust = axis_text_align_n
                                   ), axis.text.x = ggplot2::element_text(size = ggplot2::rel(size_text_relative)),
                                   legend.position = "none"
                                 ))
  if (is.null(spots) == FALSE) {
    if (is.data.frame(spots) == TRUE) {
      if (by_date == FALSE) {
        spots_date <- spots |>
          tidyr::drop_na() |>
          dplyr::mutate(
            spot_date = as.numeric(spot_date),
            activity = as.character(activity), spot_type = as.character(spot_type)
          ) |>
          dplyr::mutate(
            activity = factor(
              x = activity,
              levels = df_levels
            ), spot_date = zoo::as.Date(start_yearmon +
                                          (1 / 12) * zoo::as.yearmon(spot_date), frac = 0.5),
            end_date = as.Date(NA), wp = NA
          )
      } else {
        if (exact_date == TRUE) {
          spots_date <- spots |>
            tidyr::drop_na() |>
            dplyr::mutate(activity = as.character(activity)) |>
            dplyr::mutate(
              activity = factor(
                x = activity,
                levels = df_levels
              ), spot_date = as.Date(spot_date),
              end_date = as.Date(NA), wp = NA
            )
        } else {
          spots_date <- spots |>
            tidyr::drop_na() |>
            dplyr::mutate(activity = factor(
              x = activity,
              levels = df_levels
            ), spot_date = zoo::as.Date(zoo::as.yearmon(spot_date),
                                        frac = 0.5
            ), end_date = as.Date(NA), wp = NA)
        }
      }
      gg_gantt <- gg_gantt + ggplot2::geom_label(
        data = spots_date,
        mapping = ggplot2::aes(
          x = spot_date, y = activity,
          label = spot_type
        ), colour = "gray30", fontface = "bold",
        family = font_family, size = 3 * size_text_relative
      )
    }
  }
  if (show_vertical_lines == FALSE) {
    gg_gantt <- gg_gantt + ggplot2::theme(panel.grid.major.x = ggplot2::element_line(size = 0))
  }
  
  gg_gantt <- gg_gantt + ggplot2::theme(plot.margin = grid::unit(c(0.5, 0.75, 0.5, 0.5), "inches"))
  
  return(gg_gantt)
}
