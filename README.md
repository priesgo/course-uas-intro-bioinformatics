# course-uas-intro-bioinformatics

An introductory course to bioinformatics with an emphasis on pipeline development using workflow managers.

[![Open in Cloud Shell](https://gstatic.com/cloudssh/images/open-btn.svg)](https://console.cloud.google.com/cloudshell/editor?cloudshell_git_repo=https://github.com/priesgo/course-uas-intro-bioinformatics.git)


## Getting started

**Step 0**: create a GitHub user

**Step 1**: fork this repository

**Step 2**: setup a personal access token

GitHub security has stepped up recently, the easiest for this course would be to use personal access tokens. If you have one feel free to skip this step

Go to your GitHub settings > Developer settings > Personal Access Tokens > Tokens (classic) > Generate new token (classic)

Generate the token and make sure you copy this in a safe place for later use. This will be your GitHub password during the course.

**Step 3**: open a Google Cloud Shell (free for 50 hours weekly!) but make sure you point to your previously forked repository in the URL

If your GitHub user name is `MY_GITHUB` then you want to use this URL https://console.cloud.google.com/cloudshell/editor?cloudshell_git_repo=https://github.com/MY_GITHUB/course-uas-intro-bioinformatics.git

**Step 4**: register your user name and email in git (inside the Google Cloud Shell)

setup your git credential in the google shell:
```
git config --global user.name "FIRST_NAME LAST_NAME"
git config --global user.email "MY_NAME@example.com"
```

**NOTE**: You will need to repeat this operation every time you start the shell as it does maintain its state. Any work you need to save we will need to **commit** and **push** to GitHub
