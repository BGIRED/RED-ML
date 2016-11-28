// ============================================================================
// gzstream, C++ iostream classes wrapping the zlib compression library.
// Copyright (C) 2001  Deepak Bandyopadhyay, Lutz Kettner
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================
//
// File          : gzstream.h
// Revision      : $Revision: 1.5 $
// Revision_date : $Date: 2002/04/26 23:30:15 $
// Author(s)     : Deepak Bandyopadhyay, Lutz Kettner
// 
// Standard streambuf implementation following Nicolai Josuttis, "The 
// Standard C++ Library".
// ============================================================================

#ifndef GZSTREAM_H
#define GZSTREAM_H 1

// standard C++ with new header file names and std:: namespace
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#ifdef GZSTREAM_NAMESPACE
namespace GZSTREAM_NAMESPACE {
#endif

// ----------------------------------------------------------------------------
// Internal classes to implement gzstream. See below for user classes.
// ----------------------------------------------------------------------------

class gzstreambuf : public std::streambuf {
private:
    static const int bufferSize = 47+256;    // size of data buff
    // totals 512 bytes under g++ for igzstream at the end.

    gzFile           file;               // file handle for compressed file
    char             buffer[bufferSize]; // data buffer
    int             opened;             // open/close state of stream
    int              mode;               // I/O mode

    int flush_buffer();
public:
    gzstreambuf() : opened(0) 
	{
        setp( buffer, buffer + (bufferSize-1));
        setg( buffer + 4,     // beginning of putback area
              buffer + 4,     // read position
              buffer + 4);    // end position      
        // ASSERT: both input & output capabilities will not be used together
    }
   int is_open() 
	{ 
		return opened;
	}
    gzstreambuf* open( const char* name, int open_mode);
    gzstreambuf* close();
    ~gzstreambuf() { close(); }
    
    virtual int     overflow( int c = EOF);
    virtual int     underflow();
    virtual int     sync();
	gzstreambuf& operator=(const gzstreambuf& gzs);
};

class gzstreambase : virtual public std::ios 
{
protected:
    gzstreambuf buf;
public:
    gzstreambase() 
	{ 
		init(&buf);
	}
    gzstreambase( const char* name, int open_mode);
    ~gzstreambase();
	int is_open() 
	{ 
		return buf.is_open(); 
	}
    void open( const char* name, int open_mode);
    void close();
    gzstreambuf* rdbuf() 
	{ 
		return &buf;
	}
	gzstreambase& operator=(const gzstreambase& gzst);
};

// ----------------------------------------------------------------------------
// User classes. Use igzstream and ogzstream analogously to ifstream and
// ofstream respectively. They read and write files based on the gz* 
// function interface of the zlib. Files are compatible with gzip compression.
// ----------------------------------------------------------------------------

class igzstream : public gzstreambase, public std::istream 
{
public:
    igzstream() : std::istream( &buf) 
	{
		/*gzstreambase::init(&buf);*/
	} 
    igzstream( const char* name, int open_mode = std::ios::in)
        : gzstreambase( name, open_mode), std::istream( &buf) 
	{
		/*gzstreambase::init(&buf);*/
	}  
    gzstreambuf* rdbuf() { return gzstreambase::rdbuf(); }
    void open( const char* name, int open_mode = std::ios::in) {
        gzstreambase::open( name, open_mode);
    }
};

class ogzstream : public gzstreambase, public std::ostream
{
public:
    ogzstream() : std::ostream( &buf) 
	{
		/*gzstreambase::init(&buf);*/
	}
    ogzstream( const char* name, int mode = std::ios::out)
        : gzstreambase( name, mode), std::ostream( &buf) 
	{
		/*gzstreambase::init(&buf);*/
	}  
    gzstreambuf* rdbuf() 
	{ 
		return gzstreambase::rdbuf(); 
	}
    void open( const char* name, int open_mode = std::ios::out)
	{
        gzstreambase::open( name, open_mode);
    }
	ogzstream& operator=(const ogzstream& ogzst);
};

class ogzfstream : public std::ofstream
{
	ogzstream ogz;
public:

	ogzfstream():ogz(){}
	ogzfstream(const char* name, int mode = std::ios::out):ogz(name, mode){}
	void open(const char* name, int open_mode = std::ios::out)
	{
		ogz.open(name, open_mode);
	}
	void close()
	{
		ogz.close();
	}
	void clear()
	{
		ogz.clear();
	}
	//const ogzfstream& operator= (const ogzfstream &ogzcon);	
	
	ogzfstream & operator<<(const char * value)
	{
		ogz << value;
		return *this;
	}
	ogzfstream & operator<<(const char value)
	{
		ogz << value;
		return *this;
	}
	ogzfstream & operator<<(const double value)
	{
		ogz << value;
		return *this;
	}
	ogzfstream & operator<<(const int value)
	{
		ogz << value;
		return *this;
	}
	ogzfstream & operator<<(const float value)
	{
		ogz << value;
		return *this;
	}
	ogzfstream & operator<<(const std::string & value)
	{
		ogz << value;
		return *this;
	}
};

#ifdef GZSTREAM_NAMESPACE
} // namespace GZSTREAM_NAMESPACE
#endif

#endif // GZSTREAM_H
// ============================================================================
// EOF //

