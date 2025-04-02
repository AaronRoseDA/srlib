##### ---- ROXYGEN ---- #####
#' Create a Database Connection
#'
#' This function establishes a connection to an Azure SQL Server database
#' using the provided user credentials.
#'
#' @param user_name A string specifying the database username.
#' @param password A string specifying the database password.
#'
#' @return A database connection object (`DBI::dbConnect`) if the connection is successful.
#'         If the connection fails, returns NULL and prints an error message.
#'
#' @examples
#' # Create a database connection
#' con <- create_db_connection("your_username", "your_password")
#'
#' # Use the connection to interact with the database
#' if (!is.null(con)) {
#'   dbListTables(con)  # List available tables
#'   dbDisconnect(con)  # Close the connection when done
#' }
#'
#' @export
#####
create_db_connection <- function(user_name, password) {
  con <- tryCatch({
    dbConnect(odbc::odbc(),
              Driver = "ODBC Driver 17 for SQL Server",
              Server = "symregserver.database.windows.net",
              Database = "SymReg",
              UID = user_name,
              PWD = password,
              Encrypt = "yes",
              TrustServerCertificate = "no")
  }, error = function(e) {
    cat("Error: Unable to connect to the database. Please check your credentials.\n")
    return(NULL)
  })

  if (!is.null(con)) {
    cat("Database connection established successfully.\n")
  }

  return(con)
}

con <- dbConnect(odbc::odbc(),
                 Driver = "{ODBC Driver 18 for SQL Server}",
                 Server = "symregserver.database.windows.net",
                 Database = "SymReg",
                 Authentication = "ActiveDirectoryInteractive",
                 Encrypt = "yes",
                 TrustServerCertificate = "no",
                 Timeout = 30)





##### ---- ROXYGEN ---- #####
#' Append or Replace Data in a Database Table
#'
#' This function appends or replaces data in an existing database table.
#' If the table does not exist, it prompts the user to select from available tables.
#' The function also ensures that the `DATE_INSERTED` column is updated if required.
#'
#' @param connection A valid database connection object.
#' @param table_name A string specifying the name of the table to append or replace data in.
#' @param data A data frame containing the data to be inserted.
#' @param verbose A logical value (default: FALSE). If TRUE, prints a confirmation message upon successful insertion.
#' @param DATE_INSERTED A logical value (default: TRUE). If TRUE, updates or adds a `DATE_INSERTED` column
#'        to match the current timestamp before insertion.
#' @param mode A string specifying the insertion mode. Must be either "append" (default) or "replace".
#'        "append" adds the data to the table, while "replace" clears the table before inserting new data.
#'
#' @return NULL. The function appends or replaces data in the specified table and prints messages based on user input.
#'
#' @examples
#' # Append data to an existing table
#' df <- data.frame(ID = 6:10, Name = c("F", "G", "H", "I", "J"))
#' append_to_table(my_connection, table_name = "ExistingTable", data = df)
#'
#' # Replace all data in an existing table
#' new_data <- data.frame(ID = 1:5, Name = c("A", "B", "C", "D", "E"))
#' append_to_table(my_connection, table_name = "ExistingTable", data = new_data, mode = "replace")
#'
#' # Append data with verbose output
#' append_to_table(my_connection, table_name = "ExistingTable", data = df, verbose = TRUE)
#'
#' @export
#####
append_to_table <- function(connection, table_name, data, verbose = FALSE, DATE_INSERTED = TRUE, mode = "append") {

  # Validate mode
  if (!mode %in% c("append", "replace")) stop("Invalid mode. Use 'append' or 'replace'.")

  # Check if the table exists
  tables <- dbGetQuery(connection, "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'dbo'")
  if (!(table_name %in% tables$TABLE_NAME)) {
    cat("Table not found. Available tables:\n")
    print(tables)
    table_name <- readline("Enter the correct table name (exact match required) or press Enter to cancel: ")
    if (!(table_name %in% tables$TABLE_NAME)) {
      cat("Operation canceled.\n")
      return(NULL)
    }
  }

  # Fetch the table structure
  db_columns <- dbListFields(connection, table_name)

  # Ensure DATE_INSERTED is in both input and database schema
  if (DATE_INSERTED) {
    if ("DATE_INSERTED" %in% db_columns) {
      data$DATE_INSERTED <- Sys.time()  # Update or add DATE_INSERTED
    } else {
      data$DATE_INSERTED <- NULL  # Remove it if not in DB schema
    }
  }

  # Ensure column names match
  data <- data[, intersect(names(data), db_columns), drop = FALSE]

  # Confirm before replacing data
  if (mode == "replace") {
    confirm <- readline(sprintf("Are you sure you want to replace all data in '%s'? Type 'yes' to confirm: ", table_name))
    if (tolower(confirm) != "yes") {
      cat("Operation canceled.\n")
      return(NULL)
    }
    dbExecute(connection, paste0("DELETE FROM " , table_name))
  }

  # Append the data
  dbWriteTable(connection, table_name, data, append = TRUE, row.names = FALSE)

  if (verbose) cat("Data successfully ", mode, "ed into table: ", table_name, "\n",sep = "")
}

##### ---- ROXYGEN ---- #####
#' Create or Update a Table in a Database
#'
#' This function creates a new table in a database based on the structure of a provided data frame.
#' If the table already exists, the user is prompted to either append or replace the data.
#'
#' @param connection A valid database connection object.
#' @param table_name A string specifying the name of the table to be created or updated.
#' @param table A data frame defining the structure of the table. If the data frame contains rows, they will be inserted into the database.
#' @param DATE_INSERTED A logical value (default: TRUE). If TRUE, ensures a `DATE_INSERTED` column exists in the database
#'        and updates all its values to the current timestamp before insertion.
#'
#' @return NULL. The function creates or updates a table and prints confirmation messages.
#'
#' @examples
#' # Create a new table with data
#' df <- data.frame(ID = 1:5, Name = c("A", "B", "C", "D", "E"))
#' create_table(my_connection, table_name = "NewTable", table = df)
#'
#' # Create an empty table
#' empty_df <- data.frame(ID = integer(), Name = character())
#' create_table(my_connection, table_name = "EmptyTable", table = empty_df)
#'
#' # Append or replace data in an existing table
#' new_data <- data.frame(ID = 6:10, Name = c("F", "G", "H", "I", "J"))
#' create_table(my_connection, table_name = "NewTable", table = new_data)
#'
#' @export
#####
create_table <- function(connection, table_name, table, DATE_INSERTED = TRUE) {
  # connection <-  con
  # table_name <- "input_NORMAL"
  # table<- NORMAL_data
  # DATE_INSERTED = FALSE
  # Validate input
  if (!is.data.frame(table)) stop("The 'table' parameter must be a data frame.")

  # Check if the table exists
  tables <- dbGetQuery(connection, "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'dbo'")
  table_exists <- table_name %in% tables$TABLE_NAME

  if (table_exists) {
    # Table exists, prompt for append/replace
    confirm <- readline(sprintf("Table '%s' already exists. Append or replace? (append/replace/cancel): ", table_name))
    if (!confirm %in% c("append", "replace")) {
      cat("Operation canceled.\n")
      return(NULL)
    }
    append_to_table(connection, table_name, table, verbose = TRUE, DATE_INSERTED = DATE_INSERTED, mode = confirm)
  } else {
    # Create a new table
    # if (DATE_INSERTED && !"DATE_INSERTED" %in% names(table)) table$DATE_INSERTED <- Sys.time()
    if (nrow(table) > 0) {
      if (DATE_INSERTED && !"DATE_INSERTED" %in% names(table)) {
        table$DATE_INSERTED <- Sys.time()
      }
    } else {
      table$DATE_INSERTED <- character()
    }

    # Generate CREATE TABLE SQL
    create_stmt <- paste0("CREATE TABLE ", table_name, " (",
                          paste(sprintf("[%s] %s", names(table), ifelse(sapply(table, is.numeric), "FLOAT", "VARCHAR(255)")), collapse = ", "),
                          ");")

    # Execute table creation
    dbExecute(connection, create_stmt)

    # Insert data if available
    if (nrow(table) > 0) append_to_table(connection, table_name, table, verbose = TRUE, DATE_INSERTED = DATE_INSERTED, mode = "append")

    cat("Table '", table_name, "' created successfully.\n")
  }
}

##### ---- ROXYGEN ---- #####
#' Delete Tables or Views from a Database
#'
#' This function deletes specified tables or views from a database. If no items are provided,
#' it lists all objects of the specified type (`mode`) and prompts the user for confirmation before deletion.
#'
#' @param connection A valid database connection object.
#' @param items A character vector specifying the names of tables or views to delete.
#'        If NULL, the function retrieves all objects of the specified mode and prompts the user before deletion.
#' @param mode A string specifying the type of object to delete. Must be either "table" or "view".
#'
#' @return NULL. The function performs deletions interactively and prints confirmation messages.
#'
#' @examples
#' # Delete specific tables
#' delete_db_items(my_connection, items = c("Table1", "Table2"), mode = "table")
#'
#' # Delete specific views
#' delete_db_items(my_connection, items = c("View1", "View2"), mode = "view")
#'
#' # Delete all tables (asks for confirmation for each)
#' delete_db_items(my_connection, mode = "table")
#'
#' # Delete all views (asks for confirmation for each)
#' delete_db_items(my_connection, mode = "view")
#'
#' @export
#####
delete_db_items <- function(connection, items = NULL, mode) {

  # Validate mode
  if (!mode %in% c("table", "view")) {
    stop("Invalid mode. Use 'table' or 'view'.")
  }

  # If no items are provided, get all objects of the specified mode
  if (is.null(items) || length(items) == 0) {
    query <- if (mode == "table") {
      "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'dbo'"
    } else {
      "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.VIEWS WHERE TABLE_SCHEMA = 'dbo'"
    }
    items <- dbGetQuery(connection, query)$TABLE_NAME
  }

  # Iterate through each item
  for (item in items) {

    # Fetch first 10 rows before deletion
    preview_query <- sprintf("SELECT TOP 10 * FROM dbo.[%s]", item)
    preview_data <- tryCatch(dbGetQuery(connection, preview_query), error = function(e) NULL)

    if (!is.null(preview_data)) {
      print(preview_data)  # Show data preview
      cat(sprintf("\nAre you sure you want to delete '%s' (mode: %s)? Type 'yes' to confirm: ", item, mode))
      user_input <- readline()

      if (tolower(user_input) != "yes") {
        cat(sprintf("Skipping deletion of '%s'.\n", item))
        next  # Skip to the next item
      }
    } else {
      cat(sprintf("'%s' does not exist or is empty.\n", item))
      next
    }

    # If mode is table, drop foreign key constraints first
    if (mode == "table") {
      constraint_query <- sprintf("
        SELECT CONSTRAINT_NAME
        FROM INFORMATION_SCHEMA.TABLE_CONSTRAINTS
        WHERE CONSTRAINT_TYPE = 'FOREIGN KEY'
        AND TABLE_SCHEMA = 'dbo'
        AND TABLE_NAME = '%s'", item)

      constraints <- dbGetQuery(connection, constraint_query)

      if (nrow(constraints) > 0) {
        for (constraint in constraints$CONSTRAINT_NAME) {
          drop_constraint_query <- sprintf("ALTER TABLE dbo.[%s] DROP CONSTRAINT [%s];", item, constraint)
          try(dbExecute(connection, drop_constraint_query), silent = TRUE)
        }
      }

      # Drop the table
      delete_query <- sprintf("DROP TABLE dbo.[%s];", item)
    } else {  # Mode is 'view'
      delete_query <- sprintf("DROP VIEW dbo.[%s];", item)
    }

    # Execute deletion
    tryCatch({
      dbExecute(connection, delete_query)
      cat(sprintf("'%s' deleted successfully.\n", item))
    }, error = function(e) {
      cat(sprintf("Error deleting '%s': %s\n", item, e$message))
    })
  }

  cat("Deletion process completed.\n")
}

##### ---- ROXYGEN ---- #####
#' Retrieve a Table or Summary from a Database
#'
#' This function fetches data from a specified table in the database,
#' with options to retrieve all rows, a specified number of rows,
#' or a summary of all tables in the database.
#'
#' @param connection A valid database connection object.
#' @param table A string specifying the name of the table to retrieve.
#'        If `summary = TRUE`, this parameter is ignored.
#' @param n_rows An optional integer specifying the number of rows to retrieve.
#'        If NULL, the function returns the full table.
#' @param summary A logical value (default: FALSE). If TRUE, the function
#'        returns a list of all tables in the database, each with its
#'        first 4 rows, and ignores the `table` and `n_rows` parameters.
#' @param verbose A logical value (default: FALSE). If TRUE, prints the
#'        first few rows of the requested table or summary list for reference.
#'
#' @return A data frame containing the requested table data if `summary = FALSE`,
#'         or a named list of data frames containing the first 4 rows of each table
#'         in the database if `summary = TRUE`.
#'
#' @examples
#' # Retrieve the full table
#' df <- get_db_table(my_connection, table = "SalesData")
#'
#' # Retrieve the first 50 rows and print structure
#' df <- get_db_table(my_connection, table = "SalesData", n_rows = 50)
#'
#' # Get a summary of all tables (first 4 rows of each)
#' summary_list <- get_db_table(my_connection, summary = TRUE)
#'
#' # Retrieve table with verbose output
#' df <- get_db_table(my_connection, table = "SalesData", verbose = TRUE)
#'
#' @export
#####
get_db_table <- function(connection, table = NULL, n_rows = NULL, summary = FALSE, verbose = FALSE) {

  # If summary = TRUE, ignore table & n_rows and return all tables with first 10 rows
  if (summary) {
    tables <- dbGetQuery(connection,
                         "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'dbo'")
    table_list <- list()


    if (is.null(table))
      for (tbl in tables$TABLE_NAME) {
        query <- sprintf("SELECT TOP 4 * FROM dbo.[%s]", tbl)
        table_list[[tbl]] <- dbGetQuery(connection, query)
      } else if (table %in% tables$TABLE_NAME){
        query <- sprintf("SELECT TOP 4 * FROM dbo.[%s]", table)
        table_list[[table]] <- dbGetQuery(connection, query)
      }

    if (verbose){
      print(table_list)}
    return(table_list)
  }

  # Validate table parameter
  if (is.null(table)) {
    stop("Please specify a valid table name or set summary = TRUE.")
  }

  # Check if table exists
  tables <- dbGetQuery(connection,
                       "SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA = 'dbo'")

  if (!(table %in% tables$TABLE_NAME)) {
    cat("Invalid table name. Available tables:\n")
    print(tables)
    return(NULL)
  }

  # Construct query based on n_rows parameter
  query <- if (is.null(n_rows)) {
    sprintf("SELECT * FROM dbo.[%s]", table)
  } else {
    sprintf("SELECT TOP %d * FROM dbo.[%s]", n_rows, table)
  }

  # Fetch data
  data <- dbGetQuery(connection, query)
  if (verbose){
    print(head(data))
  }

  # If n_rows is specified, print table structure and head
  if (!is.null(n_rows)) {
    cat("\nTable Structure:\n")
    str(data)
    cat("\nTable Head:\n")
    print(head(data))
  }

  return(data)
}




create_table(con, "input_NORMAL", NORMAL_data, DATE_INSERTED = FALSE)
create_table(con, "input_HALFNORMAL", HALFNORMAL_data, DATE_INSERTED = FALSE)
create_table(con, "input_UNIFORM", UNIFORM_data, DATE_INSERTED = FALSE)


append_to_table(con, "input_NORMAL", NORMAL_data, DATE_INSERTED = FALSE)
append_to_table(con, "input_HALFNORMAL", HALFNORMAL_data, DATE_INSERTED = FALSE)
append_to_table(con, "input_UNIFORM", UNIFORM_data, DATE_INSERTED = FALSE)

get_db_table(connection = con, summary = T)
temp <- get_db_table(connection = con, table = "input_NORMAL", summary = F)



