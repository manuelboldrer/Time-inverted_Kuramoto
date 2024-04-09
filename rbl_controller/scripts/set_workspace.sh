#!/bin/bash

set -e

trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
trap 'echo "$0: \"${last_command}\" command failed with exit code $?"' ERR

# shift
OPTIND=1
NO_BUILD=false
while getopts "g:l:n" options; do
  case ${options} in
    g)
      GIT_PATH=${OPTARG}
      echo "Parsed GIT_PATH=$GIT_PATH"
      ;;
    l)
      WORKSPACE_LOCATION=${OPTARG}
      echo "Parsed WORKSPACE_LOCATION=$WORKSPACE_LOCATION"
      ;;
    n)
      NO_BUILD=true
      echo "NO_BUILD=true"
      ;;
  esac
done

[ -z "$WORKSPACE_LOCATION" ] && WORKSPACE_LOCATION="$HOME"
[ -z "$GIT_PATH" ] && GIT_PATH="$HOME"/git

# get the path to this script
APP_PATH=`dirname "$0"`
APP_PATH=`( cd "$APP_PATH" && pwd )`

#get the repository name 
REPOSITORY_PATH=`( cd "$APP_PATH/.." && pwd )`

WORKSPACE_NAME=manuel_workspace
WORKSPACE_PATH=$WORKSPACE_LOCATION/$WORKSPACE_NAME

# create the folder structure
mkdir -p $WORKSPACE_PATH/src

cd $WORKSPACE_PATH
command catkin init

echo "$0: setting up build profiles"
command catkin config --profile debug --cmake-args -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_CXX_FLAGS='-std=c++17 -Og' -DCMAKE_C_FLAGS='-Og'
command catkin profile set debug
# command catkin config --extend $MRS_WORKSPACE/devel

command catkin config --profile release --cmake-args -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_CXX_FLAGS='-std=c++17'
command catkin profile set release
# command catkin config --extend $MRS_WORKSPACE/devel

command catkin config --profile reldeb --cmake-args -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_CXX_FLAGS='-std=c++17'
command catkin profile set reldeb
# command catkin config --extend $MRS_WORKSPACE/devel

# build profile for normal installation
[ -z "$GITHUB_CI" ] && command catkin profile set reldeb

# build profile for github actions
[ ! -z "$GITHUB_CI" ] && command catkin profile set debug

echo "$0: linking ros packages to $WORKSPACE_LOCATION/$WORKSPACE"
echo "$REPOSITORY_PATH"
cd $WORKSPACE_PATH/src
ln -sf $REPOSITORY_PATH

echo "$0: building $WORKSPACE_PATH"
cd $WORKSPACE_PATH

if ! $NO_BUILD; then
  [ -z "$GITHUB_CI" ] && command catkin build
  [ ! -z "$GITHUB_CI" ] && command catkin build --limit-status-rate 0.2 --summarize
fi

num=`cat ~/.bashrc | grep "$WORKSPACE_PATH" | wc -l`
if [ "$num" -lt "1" ]; then

  # set bashrc
  echo "
source $WORKSPACE_PATH/devel/setup.bash" >> ~/.bashrc

fi
