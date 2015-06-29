/** \file hribf_buffers.h
  * 
  * \brief Handles poll2 output data files
  * 
  * The classes within this file are used to open, format,
  * and properly write data to an hribf style ldf file so
  * that it may be read by programs which read legacy ldf
  * files. Each individual buffer class (inheriting from the
  * BufferType class) is responsible for writting a different
  * type of buffer to the data file. PollOutputFile is the
  * class used to stitch all of the individual buffers together
  * into a correctly formatted output ldf file which may be
  * read by SCANOR from the hhirf upak library.
  *
  * \author Cory R. Thornsberry
  * 
  * \date May 6th, 2015
  * 
  * \version 1.1.05
*/

#ifndef HRIBF_BUFFERS_H
#define HRIBF_BUFFERS_H

#include <fstream>
#include <vector>

#define HRIBF_BUFFERS_VERSION "1.1.07"
#define HRIBF_BUFFERS_DATE "June 1st, 2015"

class BufferType{
  protected:
	int bufftype;
	int buffsize;
	int buffend;
	int zero;
	bool debug_mode;
	
	BufferType(int bufftype_, int buffsize_, int buffend_=-1);
	
	/// Returns only false if not overloaded
	virtual bool Write(std::ofstream *file_);

	/// Returns only false if not overloaded
	virtual bool Read(std::ifstream *file_);
	
  public:
	int GetBufferType(){ return bufftype; }
	
	int GetBufferSize(){ return buffsize; }
	
	int GetBufferEndFlag(){ return buffend; }
	
	bool DebugMode(){ return debug_mode; }
  
	void SetDebugMode(bool debug_=true){ debug_mode = debug_; }
};

/* The DIR buffer is written at the beginning of each .ldf file. When the file is ready
   to be closed, the data within the DIR buffer is re-written with run information. */
class DIR_buffer : public BufferType{
  private:
  	int total_buff_size;
	int run_num;
	int unknown[3];
	
  public:
	DIR_buffer();
	
	int GetTotalBufferSize(){ return total_buff_size; }
	
	int GetRunNumber(){ return run_num; }
	
	void SetRunNumber(int input_){ run_num = input_; }
	
	/* DIR buffer (1 word buffer type, 1 word buffer size, 1 word for total buffer length,
	   1 word for total number of buffers, 2 unknown words, 1 word for run number, 1 unknown word,
	   and 8186 zeros) */
	bool Write(std::ofstream *file_);

	/// Read a DIR buffer from a file. Return false if buffer has the wrong header and return true otherwise
	bool Read(std::ifstream *file_, int &number_buffers);
};

/* The HEAD buffer is written after the DIR buffer for each .ldf file. HEAD contains information
   about the run including the date/time, the title, and the run number. */
class HEAD_buffer : public BufferType{
  private:
	char facility[9];
	char format[9];
	char type[17];
	char date[17];
	char run_title[81];
	int run_num;

	void set_char_array(std::string input_, char *arr_, unsigned int size_);
	
  public:
	HEAD_buffer();
		
	char *GetFacility(){ return facility; }
	
	char *GetFormat(){ return format; }
		
	char *GetType(){ return type; }
	
	char *GetDate(){ return date; }
	
	char *GetRunTitle(){ return run_title; }
		
	int GetRunNumber(){ return run_num; }
		
	bool SetDateTime();
	
	bool SetTitle(std::string input_);
	
	void SetRunNumber(int input_){ run_num = input_; }

	/** HEAD buffer (1 word buffer type, 1 word buffer size, 2 words for facility, 2 for format, 
	  * 3 for type, 1 word separator, 4 word date, 20 word title [80 character], 1 word run number,
	  * 30 words of padding, and 8129 end of buffer words) */
	bool Write(std::ofstream *file_);

	/// Read a HEAD buffer from a file. Return false if buffer has the wrong header and return true otherwise
	bool Read(std::ifstream *file_);
};

/// The DATA buffer contains all physics data within the .ldf file
class DATA_buffer : public BufferType{
  private:
	unsigned int current_buff_pos; /// Absolute buffer position
	unsigned int buff_words_remaining; /// Absolute number of buffer words remaining
	unsigned int good_words_remaining; /// Good buffer words remaining (not counting header or footer words)

	/// DATA buffer (1 word buffer type, 1 word buffer size)
	bool open_(std::ofstream *file_);
	
  public:
	DATA_buffer(); /// 0x41544144 "DATA"

	/// Close a data buffer by padding with 0xFFFFFFFF
	bool Close(std::ofstream *file_);

	/** Get the standard data spill size for a given data file. This number is set at runtime by poll
	  * and should be the same for each and every spill in the file. Returns the spill size in words 
	  * or -1 in the event of an error. This method should be called whenever a new file is opened. */
	int GetSpillSize(std::ifstream *file_);
	
	/// Write a data spill to file
	bool Write(std::ofstream *file_, char *data_, unsigned int nWords_, int &buffs_written, int output_format_=0);
	
	/// Read a data spill from a file
	bool Read(std::ifstream *file_, char *data_, unsigned int &nWords_, unsigned int max_bytes_, bool &full_spill, bool &bad_spill, int file_format_=0);
};

/// A single EOF buffer signals the end of a run (pacman .ldf format). A double EOF signals the end of the .ldf file.
class EOF_buffer : public BufferType{	
  public:
	EOF_buffer(); /// 0x20464F45 "EOF "
	
	/// EOF buffer (1 word buffer type, 1 word buffer size, and 8192 end of buffer words)
	bool Write(std::ofstream *file_);

	/// Read an EOF buffer from a file. Return false if buffer has the wrong header and return true otherwise
	bool Read(std::ifstream *file_);
};

class PollOutputFile{
  private:
	std::ofstream output_file;
	std::string fname_prefix;
	std::string current_filename;
	std::string current_full_filename;
	DIR_buffer dirBuff;
	HEAD_buffer headBuff;
	DATA_buffer dataBuff;
	EOF_buffer eofBuff;
	int current_file_num;
	int output_format;
	int number_spills;
	bool debug_mode;
	int run_num;
	
	int current_depth;
	std::string current_directory;
	std::vector<std::string> directories;

	/// Get the formatted filename of the current file
	std::string get_filename();
	
	/// Get the full path of the current file
	bool get_full_filename(std::string &output);

	/** Overwrite the fourth word of the file with the total number of buffers and close the file
	  * Returns false if no output file is open or if the number of 4 byte words in the file is not 
	  * evenly divisible by the number of words in a buffer */
	bool overwrite_dir(int total_buffers_=-1);

	/// Initialize the output file with initial parameters
	void initialize();

  public:
	PollOutputFile();

	PollOutputFile(std::string filename_);
	
	~PollOutputFile(){ CloseFile(); }
	
	/// Get the size of the current file, in bytes.
	std::streampos GetFilesize(){ return output_file.tellp(); }
	
	/// Get the name of the current output file
	std::string GetCurrentFilename(){ return current_filename; }
	
	/// Return the total number of spills written since the current file was opened
	int GetNumberSpills(){ return number_spills; }
	
	/// Toggle debug mode
	void SetDebugMode(bool debug_=true);
	
	/// Set the output file format
	bool SetFileFormat(int format_);

	/// Set the output filename prefix
	void SetFilenamePrefix(std::string filename_);

	/// Return true if an output file is open and writable and false otherwise
	bool IsOpen(){ return (output_file.is_open() && output_file.good()); }
	
	/// Write nWords_ of data to the file
	int Write(char *data_, unsigned int nWords_);

	/** Build a data spill notification message for broadcast onto the network
	  * Return the total number of bytes in the packet upon success, and -1 otherwise */
	int BuildPacket(char *&output);

	/// Close the current file, if one is open, and open a new file for data output
	bool OpenNewFile(std::string title_, int &run_num_, std::string prefix, std::string output_dir="./", bool continueRun = false);

	std::string GetNextFileName(int &run_num_, std::string prefix, std::string output_dir, bool continueRun = false);
	
	int GetRunNumber() {return dirBuff.GetRunNumber();}

	/// Write the footer and close the file
	void CloseFile();
};

#endif
