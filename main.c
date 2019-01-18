/*
Copyright (c) 2019 Felipe Ferreira da Silva

This software is provided 'as-is', without any express or implied warranty. In
no event will the authors be held liable for any damages arising from the use of
this software.

Permission is granted to anyone to use this software for any purpose, including
commercial applications, and to alter it and redistribute it freely, subject to
the following restrictions:

  1. The origin of this software must not be misrepresented; you must not claim
     that you wrote the original software. If you use this software in a
     product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#ifdef MULTITHREAD_IMPLEMENTATION
#include <pthread.h>
#endif

#define APPLICATION_VERSION_YYYY 2019
#define APPLICATION_VERSION_MM 01
#define APPLICATION_VERSION_DD 28
#define APPLICATION_VERSION_MICRO 0

#define PATH_LENGTH_MAX 4096
#define MINIMUM_SIMILARITY_MAXIMUM 100000000
#define THREADS_MAX 64

enum state {
	STATE_NONE,
	STATE_HEADER,
	STATE_SEQUENCE
};

struct entry {
	uint32_t parent_index;
	uint32_t parent_matches;
	uint32_t parent_mismatches;
	uint32_t similarity;
	uint32_t position;
	uint32_t header_size;
	uint32_t sequence_size;
};

struct context {
	uint8_t print_help;
	uint8_t file_input_path[PATH_LENGTH_MAX];
	uint8_t file_output_path[PATH_LENGTH_MAX];
	FILE *file_input;
	FILE *file_output;
	uint32_t pattern_length;
	uint32_t minimum_similarity;
	uint8_t try_complement;
	uint32_t threads_total;
	struct entry *entries;
	uint32_t entries_total;
	uint32_t entry_query_index;
	uint8_t show_progress;
	uint32_t longest_sequence_size;
	uint32_t smallest_sequence_size;
#ifdef MULTITHREAD_IMPLEMENTATION
	pthread_t threads[THREADS_MAX];
#endif
};

#ifdef MULTITHREAD_IMPLEMENTATION
pthread_mutex_t mutex;
pthread_cond_t condition;
pthread_attr_t attributes;
#endif
#define PROGRESS_DIVISOR 1000000

static void get_sequences_distance(uint32_t pattern_length,
uint8_t *entry_query_sequence_bytes, uint32_t entry_query_sequence_size,
uint8_t *entry_target_sequence_bytes, uint32_t entry_target_sequence_size,
uint32_t *matches,
uint32_t *mismatches)
{
	uint32_t entry_query_sequence_offset;
	matches[0] = 0;
	entry_query_sequence_offset = 0;
	while (entry_query_sequence_offset + pattern_length < entry_query_sequence_size) {
		uint32_t entry_target_sequence_offset;
		entry_target_sequence_offset = 0;
		while (entry_target_sequence_offset + pattern_length < entry_target_sequence_size) {
			if (memcmp(entry_query_sequence_bytes + entry_query_sequence_offset, entry_target_sequence_bytes + entry_target_sequence_offset, pattern_length) == 0) {
				matches[0] = matches[0] + 1;
			} else {
				mismatches[0] = mismatches[0] + 1;
			}
			entry_target_sequence_offset = entry_target_sequence_offset + pattern_length;
		}
		entry_query_sequence_offset = entry_query_sequence_offset + pattern_length;
	}
}

static void process_entry_query(struct context *context, uint32_t entry_query_index, uint8_t *entry_query_sequence_bytes, uint8_t *entry_target_sequence_bytes)
{
	struct entry *entry_query;
	entry_query = &context->entries[entry_query_index];
	fseek(context->file_input, entry_query->position + entry_query->header_size + 1, SEEK_SET);
	if (fread(entry_query_sequence_bytes, 1, entry_query->sequence_size, context->file_input) == entry_query->sequence_size) {
		uint32_t entry_target_index;
		entry_target_index = entry_query_index + 1;
		while (entry_target_index < context->entries_total) {
			uint32_t matches;
			uint32_t mismatches;
			struct entry *entry_target;
			entry_target = &context->entries[entry_target_index];
			if (entry_target != entry_query) {
				fseek(context->file_input, entry_target->position + entry_target->header_size + 1, SEEK_SET);
				if (fread(entry_target_sequence_bytes, 1, entry_target->sequence_size, context->file_input) == entry_target->sequence_size) {
					get_sequences_distance(context->pattern_length,
						entry_query_sequence_bytes, entry_query->sequence_size,
						entry_target_sequence_bytes, entry_target->sequence_size,
						&matches,
						&mismatches);
					if (entry_target->parent_index == 0) {
						entry_target->parent_index = entry_query_index + 1;
						entry_target->parent_matches = matches;
						entry_target->parent_mismatches = mismatches;
					} else {
						if (matches > entry_target->parent_matches) {
							entry_target->parent_index = entry_query_index + 1;
							entry_target->parent_matches = matches;
							entry_target->parent_mismatches = mismatches;
						} else if (matches == entry_target->parent_matches && mismatches < entry_target->parent_mismatches) {
							entry_target->parent_index = entry_query_index + 1;
							entry_target->parent_matches = matches;
							entry_target->parent_mismatches = mismatches;
						}
					}
				}
			}
			entry_target_index = entry_target_index + 1;
			if (context->show_progress == 1 && (entry_target_index % (PROGRESS_DIVISOR / 1000) == 0 || entry_target_index == context->entries_total)) {
				fprintf(stdout, "\rClustering at %0.2f%% (entry %u at %06.2f%%).",
					((double)entry_query_index / (double)context->entries_total) * 100.0,
					entry_query_index,
					((double)entry_target_index / (double)context->entries_total) * 100.0);
				fflush(stdout);
			}
		}
	} else {
		puts("Failed to read sequence.");
	}
}

static void process_entries(struct context *context)
{
	uint8_t *entry_query_sequence_bytes;
	uint8_t *entry_target_sequence_bytes;
	entry_query_sequence_bytes = malloc(context->longest_sequence_size);
	entry_target_sequence_bytes = malloc(context->longest_sequence_size);
	fprintf(stdout, "\rClusterization at %0.2f%%.", 0.0);
	fflush(stdout);
	while (context->entry_query_index < context->entries_total) {
		process_entry_query(context, context->entry_query_index, entry_query_sequence_bytes, entry_target_sequence_bytes);
		context->entry_query_index = context->entry_query_index + 1;
		if (context->show_progress == 1 && context->entry_query_index == context->entries_total) {
			fprintf(stdout, "\rClusterization at %0.2f%%.",
				((double)context->entry_query_index / (double)context->entries_total) * 100.0);
			fflush(stdout);
		}
	}
	if (context->show_progress == 1) {
		fprintf(stdout, "\n");
		fflush(stdout);
	}
	free(entry_query_sequence_bytes);
	free(entry_target_sequence_bytes);
}

#ifdef MULTITHREAD_IMPLEMENTATION
static void *thread_handler(void *data)
{
	struct context *context;
	uint32_t entry_query_index;
	uint32_t entries_total;
	uint8_t *entry_query_sequence_bytes;
	uint8_t *entry_target_sequence_bytes;
	context = data;
	pthread_mutex_lock(&mutex);
	entry_query_index = context->entry_query_index;
	entries_total = context->entries_total;
	context->entry_query_index = context->entry_query_index + 1;
	entry_query_sequence_bytes = malloc(context->longest_sequence_size);
	entry_target_sequence_bytes = malloc(context->longest_sequence_size);
	pthread_mutex_unlock(&mutex);
	while (entry_query_index < entries_total) {
		process_entry_query(context, entry_query_index, entry_query_sequence_bytes, entry_target_sequence_bytes);
		pthread_mutex_lock(&mutex);
		entry_query_index = context->entry_query_index;
		context->entry_query_index = context->entry_query_index + 1;
		pthread_mutex_unlock(&mutex);
	}
	free(entry_query_sequence_bytes);
	free(entry_target_sequence_bytes);
	return NULL;
}

static void wait_threads(struct context *context)
{
	printf("Waiting threads.\n");
	while (context->entry_query_index < context->entries_total) {
		struct timespec timespec;
		timespec.tv_sec = 0;
		timespec.tv_nsec = 10000;
		nanosleep(&timespec, NULL);
		if (context->show_progress == 1 && (context->entry_query_index % (PROGRESS_DIVISOR * 100) == 0 || context->entry_query_index == context->entries_total)) {
			fprintf(stdout, "\rClusterization at %0.2f%%.", ((double)context->entry_query_index / (double)context->entries_total) * 100.0);
			fflush(stdout);
		}
	}
	if (context->show_progress == 1) {
		fprintf(stdout, "\rClusterization at %0.2f%%.\n", ((double)context->entry_query_index / (double)context->entries_total) * 100.0);
		fflush(stdout);
	}
}

static void terminate_threads(struct context *context)
{
	uint32_t i;
	i = 0;
	while (i < context->threads_total) {
		pthread_join(context->threads[i], NULL);
		i = i + 1;
	}
}

static void setup_threads(struct context *context)
{
	uint32_t i;
	printf("Starting %u threads.\n", context->threads_total);
	pthread_mutex_init(&mutex, NULL);
	pthread_cond_init(&condition, NULL);
	pthread_attr_init(&attributes);
	pthread_attr_setdetachstate(&attributes, PTHREAD_CREATE_JOINABLE);
	i = 0;
	while (i < context->threads_total) {
		if (pthread_create(&context->threads[i], &attributes, thread_handler, context) != 0) {
			printf("Failed to create thread %u.\n", i);
		}
		i = i + 1;
	}
}
#endif

static void read_file_input_positions(struct context *context)
{
	uint8_t state;
	uint32_t file_position;
	uint32_t file_size;
	uint32_t progress_divisor;
	uint8_t byte;
	uint32_t sequence_size;
	struct entry *entry;
	uint32_t entry_index;
	state = 0;
	sequence_size = 0;
	entry = NULL;
	entry_index = 0;
	fseek(context->file_input, 0, SEEK_END);
	file_size = ftell(context->file_input);
	progress_divisor = PROGRESS_DIVISOR;
	fseek(context->file_input, 0, SEEK_SET);
	while (fread(&byte, 1, 1, context->file_input) == 1) {
		if (state == STATE_NONE) {
			if (byte == '>') {
				state = STATE_HEADER;
				entry = &context->entries[entry_index];
				entry->position = ftell(context->file_input);
				entry_index = entry_index + 1;
			}
		} else if (state == STATE_HEADER) {
			if (byte == '\n') {
				state = STATE_SEQUENCE;
				entry->header_size = ftell(context->file_input) - entry->position - 1;
			}
		} else if (state == STATE_SEQUENCE) {
			if (byte == '>') {
				state = STATE_HEADER;
				entry->sequence_size = sequence_size;
				entry = &context->entries[entry_index];
				entry->position = ftell(context->file_input);
				entry_index = entry_index + 1;
				sequence_size = 0;
			} else {
				if (byte != '\n') {
					sequence_size = sequence_size + 1;
				}
			}
		}
		file_position = ftell(context->file_input);
		if (context->show_progress == 1 && (file_position % progress_divisor == 0 || file_position >= file_size)) {
			fprintf(stdout, "\rChecking the position of the entries in the file input at %0.2f%%.", ((double)file_position / (double)file_size) * 100.0);
			fflush(stdout);
		}
	}
	if (context->show_progress == 1) {
		fprintf(stdout, "\n");
		fflush(stdout);
	}
}

static void read_file_input_entries(struct context *context)
{
	uint8_t state;
	uint32_t file_position;
	uint32_t file_size;
	uint32_t progress_divisor;
	uint8_t byte;
	uint32_t sequence_size;
	state = 0;
	sequence_size = 0;
	fseek(context->file_input, 0, SEEK_END);
	file_size = ftell(context->file_input);
	progress_divisor = PROGRESS_DIVISOR;
	fseek(context->file_input, 0, SEEK_SET);
	while (fread(&byte, 1, 1, context->file_input) == 1) {
		if (state == STATE_NONE) {
			if (byte == '>') {
				state = STATE_HEADER;
				context->entries_total = context->entries_total + 1;
			}
		} else if (state == STATE_HEADER) {
			if (byte == '\n') {
				state = STATE_SEQUENCE;
			}
		} else if (state == STATE_SEQUENCE) {
			if (byte == '>') {
				if (sequence_size > context->longest_sequence_size) {
					context->longest_sequence_size = sequence_size;
				}
				sequence_size = 0;
				context->entries_total = context->entries_total + 1;
				state = STATE_HEADER;
			} else {
				if (byte != '\n') {
					sequence_size = sequence_size + 1;
				}
			}
		}
		file_position = ftell(context->file_input);
		if (context->show_progress == 1 && (file_position % progress_divisor == 0 || file_position >= file_size)) {
			fprintf(stdout, "\rChecking the entries in the file input at %0.2f%%.", ((double)file_position / (double)file_size) * 100.0);
			fflush(stdout);
		}
	}
	if (context->show_progress == 1) {
		fprintf(stdout, "\n");
		fflush(stdout);
	}
	printf("Total of %u sequences.\n", context->entries_total);
	printf("Longest sequence: %u characters.\n", context->longest_sequence_size);
	printf("Allocating %lu bytes.\n", sizeof(struct entry) * context->entries_total);
	context->entries_total = context->entries_total;
	context->entries = malloc(sizeof(struct entry) * context->entries_total);
	memset(context->entries, 0, sizeof(struct entry) * context->entries_total);
}

static void print_help(void)
{
	printf("SC (Sequence Clusterizer) version %04u.%02u.%02u.%u.\n",
		APPLICATION_VERSION_YYYY,
		APPLICATION_VERSION_MM,
		APPLICATION_VERSION_DD,
		APPLICATION_VERSION_MICRO);
	puts("Created by Felipe Ferreira da Silva.");
	puts("");
	puts("Usage:");
	puts("  sc [1] [2] [3] [4] [5] [6] [7]");
	puts("");
	puts("  Argument 1: Path to input file.");
	puts("  Argument 2: Path to output file.");
	puts("  Argument 3: Pattern length (larger than 0 and smaller than the smallest sequence).");
	puts("  Argument 4: Minimum sequence similarity (0-100000000).");
	puts("  Argument 5: Try complement sequence (0 for no, 1 for yes).");
	puts("  Argument 6: Total of dedicated threads to use (0 for no dedicated thread, maximum of 64).");
	puts("  Argument 7: Show progress (0 for no, 1 for yes).");
	puts("");
	puts("File input must be in single-line FASTA format.");
}

int main(int argc, char **args)
{
	uint8_t success;
	int32_t status;
	struct context *context;
	success = 1;
	context = malloc(sizeof(struct context));
	if (context == NULL) {
		success = 0;
		goto done;
	}
	memset(context, 0, sizeof(struct context));
	/* Arguments */
	if (argc == 8) {
		strncat((char *)context->file_input_path, args[1], PATH_LENGTH_MAX - 1);
		strncat((char *)context->file_output_path, args[2], PATH_LENGTH_MAX - 1);
		context->pattern_length = atoi(args[3]);
		context->minimum_similarity = atoi(args[4]);
		context->try_complement = atoi(args[5]);
		context->threads_total = atoi(args[6]);
		context->show_progress = atoi(args[7]);
	} else {
		context->print_help = 1;
	}
	/* Process */
	if (context->print_help == 1) {
		print_help();
		goto freecontext;
	}
	printf("  Input: \"%s\"\n", context->file_input_path);
	printf("  Output: \"%s\"\n", context->file_output_path);
	printf("  Pattern length: %u\n", context->pattern_length);
	printf("  Minimum sequence similarity: %u (%0.2f%%)\n",
		context->minimum_similarity,
		((double)context->minimum_similarity / (double)MINIMUM_SIMILARITY_MAXIMUM) * 100.0);
	printf("  Try complement sequence: %u\n", context->try_complement);
	printf("  Total of dedicated threads: %u\n", context->threads_total);
	context->file_input = fopen((char *)context->file_input_path, "r");
	if (context->file_input == NULL) {
		puts("Failed to open the input file for reading.");
		success = 0;
		goto freecontext;
	}
	context->file_output = fopen((char *)context->file_output_path, "w");
	if (context->file_output == NULL) {
		puts("Failed to open the output file for writing.");
		success = 0;
		goto freecontext;
	}
	if (context->pattern_length < 1) {
		puts("The minimum length of the pattern is 1.");
		success = 0;
		goto freecontext;
	}
	if (context->minimum_similarity >= MINIMUM_SIMILARITY_MAXIMUM) {
		puts("The minimum similarity is larger than allowed.");
		success = 0;
		goto freecontext;
	}
	if (context->try_complement > 1) {
		puts("The argument to try complement sequence must be 0 or 1.");
		success = 0;
		goto freecontext;
	}
	if (context->threads_total > 64) {
		puts("The total of dedicated threads is larger than allowed.");
		success = 0;
		goto freecontext;
	}
	read_file_input_entries(context);
	read_file_input_positions(context);
#ifdef MULTITHREAD_IMPLEMENTATION
	if (context->threads_total > 0) {
		setup_threads(context);
		wait_threads(context);
		terminate_threads(context);
	} else {
		process_entries(context);
	}
#else
	process_entries(context);
#endif
	/* Write output */
freecontext:
	if (context->entries != NULL) {
		free(context->entries);
	}
	if (context->file_input != NULL) {
		fclose(context->file_input);
	}
	if (context->file_output != NULL) {
		fclose(context->file_output);
	}
	free(context);
done:
	if (success == 1) {
		status = EXIT_SUCCESS;
	} else {
		status = EXIT_FAILURE;
	}
	return status;
}
