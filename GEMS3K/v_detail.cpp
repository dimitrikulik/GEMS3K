#include <cstring>
#include <istream>
#include "v_detail.h"

void strip(std::string& str)
{
  std::string::size_type pos1 = str.find_first_not_of(' ');
  std::string::size_type pos2 = str.find_last_not_of(' ');
  str = str.substr(pos1 == std::string::npos ? 0 : pos1,
    pos2 == std::string::npos ? str.length() - 1 : pos2 - pos1 + 1);
}

void replace( std::string& str, const char* old_part, const char* new_part)
{
    size_t pos = str.find( old_part ); //rfind( old_part );
    if( pos != std::string::npos )
    {
        std::string res(str.substr(0, pos));
        res += new_part;
        res += str.substr( pos+strlen(old_part));
        str = res;
    }
}

void replaceall(std::string& str, char ch1, char ch2)
{
  for(size_t ii=0; ii<str.length(); ii++ )
   if( str[ii] == ch1 )
            str[ii] = ch2;
}

/// read string as: "<characters>"
std::istream& f_getline(std::istream& is, std::string& str, char delim)
{
    char ch;
    is.get(ch);
    str="";

    while( is.good() && ( ch==' ' || ch=='\n' || ch== '\t') )
        is.get(ch);
    if(ch == '\"')
        is.get(ch);
    while( is.good() &&  ch!=delim && ch!= '\"' )
    {
        str += ch;
        is.get(ch);
    }
    while( is.good() &&  ch!=delim )
            is.get(ch);

   return is;
}

std::string
u_makepath(const std::string& dir,
           const std::string& name, const std::string& ext)
{
    std::string Path(dir);
    if( dir != "")
      Path += "/";
    Path += name;
    Path += ".";
    Path += ext;

    return Path;
}


void
u_splitpath(const std::string& Path, std::string& dir,
            std::string& name, std::string& ext)
{
    size_t pos = Path.rfind("/");
    if( pos != std::string::npos )
    {
        dir = Path.substr(0, pos);
        pos++;
    }
    else
    {
        dir = "";
        pos = 0;
    }
    size_t pose = Path.rfind(".");
    if( pose != std::string::npos )
    {
        ext = Path.substr( pose+1, std::string::npos );
        name = Path.substr(pos, pose-pos);
    }
    else
    {
        ext = "";
        name = Path.substr(pos, std::string::npos);
    }
}
