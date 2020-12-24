#!/bin/bash -e

# Install brew if required.
if [[ $(command -v brew) == "" ]]; then
  echo "Installing Hombrew"
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
fi

brew install capnp \
             coreutils \
             eigen \
             ffmpeg \
             glfw \
             libarchive \
             libusb \
             libtool \
             llvm \
             pyenv \
             qt5 \
             zeromq

if [[ $SHELL == "/bin/zsh" ]]; then
  RC_FILE="$HOME/.zshrc"
elif [[ $SHELL == "/bin/bash" ]]; then
  RC_FILE="$HOME/.bash_profile"
fi

if [ -z "$OPENPILOT_ENV" ] && [ -n "$RC_FILE" ] && [ -z "$CI" ]; then
  OP_DIR=$(git rev-parse --show-toplevel)
  echo "source $OP_DIR/tools/openpilot_env.sh" >> $RC_FILE
  source $RC_FILE
  echo "Added openpilot_env to RC file: $RC_FILE"
fi

pyenv install -s 3.8.2
pyenv global 3.8.2
pyenv rehash
eval "$(pyenv init -)"

pip install pipenv==2020.8.13
pipenv install --system --deploy
