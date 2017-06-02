#!/bin/bash
##########################################################

# F_Time
# Stores data_time to dateTime global variable
function F_Time {
  dateTime=$(date +"%Y-%m-%d_%H-%M")
}


F_Time
message="Update "$dateTime


git add .
git commit -m "$message"
git push -u origin master




##########################################################
