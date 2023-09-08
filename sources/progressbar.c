/*
 * progressbar.c
 *
 *  Created on: 25/11/2016
 *      Author: Joan Jené
 */

#include "progressbar.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/**
 * This function shows a progress bar like this:
 *
 *      Progress: [..............................100%]
 *
 * @param current is the current iteration number
 * @param total is the last iteration number
 * @param bar_size is the length of the progress bar.
 *
 * Usage:
 * Call this function in a loop like this:
 *
 *      for (i = 0; i < 500; i++) {
 *          ShowProgressBar(i, 500, 30);
 *      }
 *      printf("\n");
 */
void ShowProgressBar(double current, double total, int bar_size) {
	double percent   = current / (total-1) * 100;					/* This is the percentage value. */
	int percent_len  = (percent<9)?1:((percent<99)?2:2);    		/* This is the number of digits of the percentage value. */
	int left         = (int)percent * bar_size / 100;       		/* This is the left space of the percentage filled with dots. */
	int right        = bar_size + 2 - left - percent_len;   		/* This is the right space of the percentage filled with spaces */
	                                                        		/* +2 means that it starts with one digit and it ends with three digits. */

    char *dots = (char *)malloc(sizeof(char) * (left + 1));         /* allocate space for the left side of the percentage value. */
    if (dots != NULL) {
		char *spaces = (char *)malloc(sizeof(char) * (right + 1));  /* allocate space for the right side of the percentage value. */
		if (spaces != NULL) {
			memset(dots, '.', left);                        		/* fill the left side of the percentage with dots. */
			dots[left] = '\x0';
			memset(spaces, ' ', right);                     		/* fill the right side of the percentage with spaces. */
			spaces[right] = '\x0';

			printf("\r\tProgress: [%s%d%%%s]",              		/* draw the progress bar. */
				   dots, (int)percent, spaces);

			fflush(stdout);                                 		/* force to print immediately */

			free(spaces);								    		/* free the allocated memory */
			spaces = NULL;
		}
		free(dots);                                         		/* free the allocated memory */
		dots = NULL;
    }
}

