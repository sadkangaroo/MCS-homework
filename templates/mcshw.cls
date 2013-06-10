\NeedsTeXFormat{LaTeX2e}[2009/01/01]
\ProvidesClass{mcshw}[2012/01/01 v0.1 mcshw for Mathematics in Computer Science homework]

\DeclareOption*{
	\PassOptionsToClass{\CurrentOption}{article}
}
\ProcessOptions\relax
\LoadClass{article}

\RequirePackage{mcshw_fonts}
\RequirePackage{mcshw_title}
\RequirePackage{amsmath, amsbsy, amssymb}

\DeclareMathOperator{\lcm}{lcm}

\setlength\oddsidemargin{-0.05in}
\setlength\evensidemargin{-0.05in}
\setlength\textheight{8.0in} 
\setlength\topmargin{-0.1in}
\setlength\textwidth{6.3in}
\addtolength\footskip{0.3in}

\newenvironment{solution} {

	\vspace{2mm}
	{\bfseries\slshape Solution:~}
	\vspace{2mm}
	\item[]
} {
	\vspace{5mm}

}

\newenvironment{proof} {
	
	\vspace{2mm}
	{\bfseries\slshape Proof:~}
	\vspace{2mm}
	\begin{list}{}{\setlength\listparindent{2em} \setlength\itemindent{2em}}
	\item[]
} {
	\end{list}
	\begin{flushright}
	$\square$
	\end{flushright}
	\vspace{5mm}

}

\newenvironment{construction} {

	\vspace{2mm}
	{\bfseries\slshape Construction:~}
	\vspace{2mm}
	\begin{list}{}{\setlength\listparindent{2em} \setlength\itemindent{2em}}
	\item[]
} {
	\end{list}
	\vspace{5mm}

}


