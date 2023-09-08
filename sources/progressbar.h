/*
 * progressbar.h
 *
 *  Created on: 25/11/2016
 *      Author: Joan Jené
 */

#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

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
	void ShowProgressBar(double current, double total, int bar_size);

#endif /* PROGRESSBAR_H */
