"boa.license" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   ver <- boa.version()
   cat("\nBayesian Output Analysis Program (BOA)\n",
       "Version ", paste(c(ver$major, ver$minor, ver$revision), collapse = "."),
       " for ", ver$system, "\n",
       "Copyright (c) 2007 Brian J. Smith <brian-j-smith@uiowa.edu>\n\n",
       "This program is free software; you can redistribute it and/or\n",
       "modify it under the terms of the GNU General Public License\n",
       "as published by the Free Software Foundation; either version 2\n",
       "of the License or any later version.\n\n",
       "This program is distributed in the hope that it will be useful,\n",
       "but WITHOUT ANY WARRANTY; without even the implied warranty of\n",
       "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n",
       "GNU General Public License for more details.\n\n",
       "For a copy of the GNU General Public License write to the Free\n",
       "Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,\n",
       "MA  02111-1307, USA, or visit their web site at\n",
       "http://www.gnu.org/copyleft/gpl.html\n\n", sep = "")
   invisible()
}
