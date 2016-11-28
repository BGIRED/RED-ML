/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-7-29
 *CLASS NAME: SamCtrl
 *FUNCTION : a class that can control sam/bam format file. It is convenient to control sam/bam 
 *			format file like open, read, write, close and so on.
 *FILE NAME : SamCtrl.h
 *UPDATE DATE : 2010-7-30
 *UPDATE BY : Bill Tang
 *******************************************************************************
 */
#ifndef SAMCTRL_H
#define SAMCTRL_H

#pragma once

#include <cmath>
#include <string>
#include "sam/bam.h"
#include "sam/faidx.h"
#include "sam/knetfile.h"
#include "sam/sam.h"
#include "sam/sam_view.h"
#include "sam/zlib.h"
#include "sam/glf.h"
#include "sam/kstring.h"
#include "sam/razf.h"
#include "sam/sam_header.h"
#include "sam/zconf.h"

class SamCtrl
{
private:
	std::string m_in_path; // the sam/bam in file's path.
	std::string m_out_path; // the sam/bam out file's path.
	samfile_t *m_in; // sam/bam file handler, use as reading.
	samfile_t *m_out; // sam/bam file handler, use as writing.
	std::string m_in_mode; // reading mode
	std::string m_out_mode; // wrinting mode
	char *m_fn_list; // sam/bam file list.
	char *m_s; // an tmp varible
	bam1_t *m_b; // a bam constructor.
public:
	SamCtrl();
	~SamCtrl(void);

	bool open(const char *path, const char *mode);
	void close();
	int readline(std::string &line);
	bool isOpened();
};

#endif /* SAMCTRL_H_ */
