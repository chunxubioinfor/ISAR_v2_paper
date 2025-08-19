#!/bin/bash
# Usage: ./connect_to_github.sh <github_username> <repo_name>

USERNAME=$1
REPO=$2

if [ -z "$USERNAME" ] || [ -z "$REPO" ]; then
echo "Usage: ./connect_to_github.sh <github_username> <repo_name>"
exit 1
fi

# Initialize git repo if not already
git init

# Add remote (replace if it already exists)
git remote remove origin 2>/dev/null
git remote add origin https://github.com/$USERNAME/$REPO.git

# Stage all files
git add .

# Commit
git commit -m "Initial commit"

# Push to master (default branch)
git branch -M master
git push -u origin master

echo "✅ Project connected to GitHub at https://github.com/$USERNAME/$REPO"