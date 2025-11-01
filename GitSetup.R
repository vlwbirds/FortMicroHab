install.packages("usethis")
system("git --version")
system('git config --global user.name "vlwbirds"')
system('git config --global user.email "vlwbirds@gmail.com"')
system("git config --list")
library(gitcreds)
gitcreds::gitcreds_set()
library(usethis)
usethis::use_git_config(user.name = "vlwbirds", user.email = "vlwbirds@gmail.com")
usethis::use_git()
usethis::use_github()

# 1) Make sure Git exists on PATH (just to be safe)
system("git --version")

# 2) Update the low-level pieces first
install.packages(c("curl", "openssl"))

# 3) Update the HTTP/auth stack usethis relies on
install.packages(c("httr2", "gh", "gert", "gitcreds", "usethis"))
install.packages("httr2")

# If you use pak (faster dependency resolution), you can instead do:
# install.packages("pak"); pak::pak(c("curl","openssl","httr2","gh","gert","gitcreds","usethis"))

# 4) Restart R *completely* (Session ➜ Restart R)

# 5) Quick diagnostics
usethis::git_sitrep()

# 6) If you haven’t stored your PAT yet:
# usethis::create_github_token()  # opens browser to make a token
# gitcreds::gitcreds_set()        # paste token once

# 7) Now connect this project to GitHub
library(httr2)
usethis::use_github(protocol = "https", private = FALSE)

lib <- .libPaths()[1]
list.files(lib, pattern = "^00LOCK", full.names = TRUE)        # see locks
unlink(file.path(lib, "00LOCK-httr2"), recursive = TRUE, force = TRUE)
# also remove any half-installed httr2 dir:
unlink(file.path(lib, "httr2"), recursive = TRUE, force = TRUE)
install.packages(c("curl","openssl"), repos = "https://cloud.r-project.org", type = "binary")
install.packages("httr2", dependencies = TRUE, repos = "https://cloud.r-project.org", type = "binary")

library(httr2); packageVersion("httr2")
install.packages(c("gh","gert","gitcreds","usethis"), repos = "https://cloud.r-project.org", type = "binary")

# restart R once more (Session → Restart R), then:
usethis::git_sitrep()
# If needed: usethis::create_github_token(); gitcreds::gitcreds_set()

usethis::use_github(protocol = "https", private = FALSE)

