118 export PATH=/home/chalabaj/skripty/launchers:$PATH
119
120 bind '"\e[A"':history-search-backward
121 bind '"\e[B"':history-search-forward
122
123 # some more ls aliases
124 alias d='ls -lrt'
125 alias q='qstat -u chalabaj'
126 alias w='qstat -f | head -n290'
127 alias e='qstat -g c'
128 alias xmg='xmgrace'
129 alias a='ls -A'
130 alias s='ls -CF'
131 alias l='ls -lrt'
132 #barva
133 PS1='${debian_chroot:+($debian_chroot)}\[\033[00;37m\]\u@\h\[\033[01;00m\]:\[\033[01;33m\]\w\[\033[00m\]/ $ '
    force_color_prompt=yes
 15 # append to the history file, don't overwrite it
 16 shopt -s histappend
 17
 18 # for setting history length see HISTSIZE and HISTFILESIZE in bash(1)
 19 HISTSIZE=100000
 20 HISTFILESIZE=200000
 21
 22 # check the window size after each command and, if necessary,
 23 # update the values of LINES and COLUMNS.
 24 shopt -s checkwinsize
87 # some more ls aliases
 88 alias ll='ls -l'
 89 alias la='ls -A'
 90 alias l='ls -CF'

 23 # update the values of LINES and COLUMNS.
 24 shopt -s checkwinsize
# Alias's for archives
alias mktar='tar -cvf'
alias mkbz2='tar -cvjf'
alias mkgz='tar -cvzf'
alias untar='tar -xvf'
alias unbz2='tar -xvjf'
alias ungz='tar -xvzf'



extract () {
	for archive in $*; do
		if [ -f $archive ] ; then
			case $archive in
				*.tar.bz2)   tar xvjf $archive    ;;
				*.tar.gz)    tar xvzf $archive    ;;
				*.bz2)       bunzip2 $archive     ;;
				*.rar)       rar x $archive       ;;
				*.gz)        gunzip $archive      ;;
				*.tar)       tar xvf $archive     ;;
				*.tbz2)      tar xvjf $archive    ;;
				*.tgz)       tar xvzf $archive    ;;
				*.zip)       unzip $archive       ;;
				*.Z)         uncompress $archive  ;;
				*.7z)        7z x $archive        ;;
				*)           echo "don't know how to extract '$archive'..." ;;
			esac
		else
			echo "'$archive' is not a valid file!"
		fi
	done
}
