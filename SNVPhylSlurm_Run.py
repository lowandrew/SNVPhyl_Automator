import SNVPhylSlurm

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--force", action="store_true",
                        help="Don't ask to update redmine api key")

    args = parser.parse_args()
    automator = SNVPhylSlurm.Automate(args.force)

    # try to run the automation tasks, if an error occurs print it
    try:
        automator.timed_retrieve()
    except Exception as e:
        import traceback
        automator.timelog.time_print("[Error] Dumping...\n%s" % traceback.format_exc())
        raise
